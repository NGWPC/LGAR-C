#ifndef BMI_LGAR_CXX_INCLUDED
#define BMI_LGAR_CXX_INCLUDED


#include <stdio.h>
#include <string>
#include <cstring>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <iostream>
#include <unordered_map>
#include "../bmi/bmi.hxx"
#include "../include/bmi_lgar.hxx"
#include "../include/all.hxx"
#include "../include/Logger.hpp"

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>

std::stringstream bmilgar_ss("");
static int lasam_state_validation_log_count = 0;
static const int lasam_state_validation_log_limit = 50;
static int GLOBAL_ERROR_COUNT = 0;
static const int GLOBAL_ERROR_LIMIT = 50;

static std::unordered_map<model_state*, int> g_lasam_export_capacity;

static inline void validate_wetting_front_state(
    const char* stage,
    struct model_state* state
) {
  if (state == nullptr || state->head == nullptr) {
    return;
  }

  // ---- LOG CONTROLS ----
  static int global_log_count = 0;
  static const int global_log_limit = 2000;   // HARD CAP
  static const int log_every_n_timesteps = 50; // sampling
  static const bool only_log_bad = true;       // key switch
  // ----------------------

  // skip logging unless:
  if (global_log_count >= global_log_limit) {
    return;
  }

  if (state->lgar_bmi_params.timesteps % log_every_n_timesteps != 0) {
    return;
  }

  // optional: only log certain stages
  if (strcmp(stage, "before_bmi_export") != 0 &&
      strcmp(stage, "after_insert_water") != 0) {
    return;
  }

  struct wetting_front* current = state->head;
  int idx = 0;
  double prev_depth = -1.0;

  while (current != NULL) {

    int layer_num = current->layer_num;
    int soil_num = -1;
    double theta_r = -9999.0;
    double theta_e = -9999.0;

    if (layer_num >= 1 && layer_num <= state->lgar_bmi_params.num_layers) {
      soil_num = state->lgar_bmi_params.layer_soil_type[layer_num];
      if (soil_num >= 1 && soil_num <= state->lgar_bmi_params.num_soil_types) {
        theta_r = state->soil_properties[soil_num].theta_r;
        theta_e = state->soil_properties[soil_num].theta_e;
      }
    }

    bool bad = false;

    if (!std::isfinite(current->theta) ||
        !std::isfinite(current->psi_cm) ||
        !std::isfinite(current->depth_cm) ||
        !std::isfinite(current->dzdt_cm_per_h) ||
        !std::isfinite(current->K_cm_per_h)) {
      bad = true;
    }

    if (current->depth_cm < 0.0) bad = true;
    if (prev_depth > current->depth_cm) bad = true;

    if (theta_e > 0.0) {
      if (current->theta < theta_r - 1e-10 || current->theta > theta_e + 1e-10) {
        bad = true;
      }
    } else {
      if (current->theta < 0.0 || current->theta > 1.0) {
        bad = true;
      }
    }

    // ---- LOGGING ----
    if (!only_log_bad || bad) {

      std::stringstream msg;
      msg << "LASAM WF"
          << " stage=" << stage
          << " idx=" << idx
          << " theta=" << current->theta
          << " depth_cm=" << current->depth_cm
          << " psi_cm=" << current->psi_cm
          << " dzdt=" << current->dzdt_cm_per_h
          << " K=" << current->K_cm_per_h
          << " bad=" << bad
          << " t=" << state->lgar_bmi_params.timesteps;

      std::string s = msg.str();

      std::cerr << s << std::endl;
      LOG(s, LogLevel::INFO);

      global_log_count++;

      if (global_log_count >= global_log_limit) {
        std::cerr << "LASAM LOG LIMIT REACHED\n";
        break;
      }
    }

    if (bad) {
      return; // stop at first failure
    }

    prev_depth = current->depth_cm;
    current = current->next;
    idx++;
  }
}

// default verbosity is set to 'none' other option 'high' or 'low' needs to be specified in the config file
string verbosity="none";

/**
 * @brief Delete dynamic arrays allocated in Initialize() and held by this object
 * 
 */
BmiLGAR::~BmiLGAR(){
  delete [] giuh_ordinates;
  delete [] giuh_runoff_queue;
}

/* The `head` pointer stores the address in memory of the first member of the linked list containing
   all the wetting fronts. The contents of struct wetting_front are defined in "all.h" */

void BmiLGAR::
Initialize (std::string config_file)
{
    // Initialize the Error, Warning and Trapping System
#ifdef EWTS_HAVE_NGEN_BRIDGE    
  EwtsInit(EWTS_ID_LASAM, true);
#else
  EwtsInit(EWTS_ID_LASAM, false);
#endif

  LOG("Inside BmiLGAR::Initialize \n", LogLevel::INFO);  
  if (config_file.compare("") != 0 ) {
    this->state = new model_state;
    state->head = NULL;
    state->state_previous = NULL;
    lgar_initialize(config_file, state);
  }

  num_giuh_ordinates = state->lgar_bmi_params.num_giuh_ordinates;

  /* giuh ordinates are static and read in the lgar.cxx, and we need to have a copy of it to pass to
     giuh.cxx, so allocating/copying here*/

  giuh_ordinates = new double[num_giuh_ordinates];
  giuh_runoff_queue = new double[num_giuh_ordinates+1];

  for (int i=0; i<num_giuh_ordinates;i++){
    giuh_ordinates[i] = state->lgar_bmi_params.giuh_ordinates[i+1]; // note lgar uses 1-indexing
  }

  for (int i=0; i<=num_giuh_ordinates;i++){
    giuh_runoff_queue[i] = 0.0;
  }

}

/**
 * @brief Allocate (or reallocate) storage for soil parameters
 * 
 */
void BmiLGAR::realloc_soil(){

  if (state == nullptr) {
    return;
  }

  int needed = state->lgar_bmi_params.num_wetting_fronts;
  if (needed <= 0) {
    needed = 1;
  }

  int &capacity = g_lasam_export_capacity[state];

  // Keep the same pointer if the current capacity is enough.
  // This avoids invalidating any downstream cached pointer.
  if (capacity >= needed &&
      state->lgar_bmi_params.soil_depth_wetting_fronts != nullptr &&
      state->lgar_bmi_params.soil_moisture_wetting_fronts != nullptr) {
    return;
  }

  // Grow only; never shrink during the run.
  int new_capacity = needed;
  if (capacity > 0) {
    new_capacity = std::max(needed, capacity * 2);
  }

  double *new_depth = new double[new_capacity]();
  double *new_moisture = new double[new_capacity]();

  if (state->lgar_bmi_params.soil_depth_wetting_fronts != nullptr &&
      state->lgar_bmi_params.soil_moisture_wetting_fronts != nullptr &&
      capacity > 0) {

    int copy_count = std::min(capacity, new_capacity);

    for (int i = 0; i < copy_count; i++) {
      new_depth[i] = state->lgar_bmi_params.soil_depth_wetting_fronts[i];
      new_moisture[i] = state->lgar_bmi_params.soil_moisture_wetting_fronts[i];
    }

    delete [] state->lgar_bmi_params.soil_depth_wetting_fronts;
    delete [] state->lgar_bmi_params.soil_moisture_wetting_fronts;
  }

  state->lgar_bmi_params.soil_depth_wetting_fronts = new_depth;
  state->lgar_bmi_params.soil_moisture_wetting_fronts = new_moisture;
  capacity = new_capacity;
}

/*
  This is the main function calling lgar subroutines for creating, moving, and merging wetting fronts.
  Calls to AET and mass balance module are also happening here
  If the model's timestep is smaller than the forcing's timestep then we take subtimesteps inside the subcycling loop
*/

void BmiLGAR::
Update()
{
  if (verbosity.compare("none") != 0) {
    bmilgar_ss <<"---------------------------------------------------------\n";
    bmilgar_ss <<"|****************** LASAM BMI Update... ******************|\n";
    bmilgar_ss <<"---------------------------------------------------------\n";
    LOG(bmilgar_ss.str(), LogLevel::INFO); bmilgar_ss.str("");
  }

  static int GLOBAL_ERROR_COUNT = 0;
  static const int GLOBAL_ERROR_LIMIT = 50;
  static bool FIRST_BAD_FOUND = false;

  double mm_to_cm = 0.1; // unit conversion
  double mm_to_m = 0.001;

  if (state->lgar_bmi_params.is_invalid_soil_type) {
    // add to mass balance accumulated variables
    state->lgar_mass_balance.volprecip_cm  += state->lgar_bmi_input_params->precipitation_mm_per_h * mm_to_cm;
    state->lgar_mass_balance.volin_cm       = 0.0;
    state->lgar_mass_balance.volon_cm       = 0.0;
    state->lgar_mass_balance.volend_cm      = state->lgar_mass_balance.volstart_cm;
    state->lgar_mass_balance.volAET_cm      = 0.0;
    state->lgar_mass_balance.volrech_cm     = 0.0;
    state->lgar_mass_balance.volrunoff_cm  += state->lgar_bmi_input_params->precipitation_mm_per_h * mm_to_cm;
    state->lgar_mass_balance.volQ_cm       += state->lgar_bmi_input_params->precipitation_mm_per_h * mm_to_cm;
    state->lgar_mass_balance.volQ_gw_cm     = 0.0;
    state->lgar_mass_balance.volPET_cm      = 0.0;
    state->lgar_mass_balance.volrunoff_giuh_cm  = 0.0;
    state->lgar_mass_balance.volchange_calib_cm = 0.0;

    // converted values, a struct local to the BMI and has bmi output variables
    bmi_unit_conv.mass_balance_m        = 0.0;
    bmi_unit_conv.volprecip_timestep_m  = state->lgar_bmi_input_params->precipitation_mm_per_h * mm_to_m;
    bmi_unit_conv.volin_timestep_m      = 0.0;
    bmi_unit_conv.volend_timestep_m     = 0.0;
    bmi_unit_conv.volAET_timestep_m     = 0.0;
    bmi_unit_conv.volrech_timestep_m    = 0.0;
    bmi_unit_conv.volrunoff_timestep_m  = state->lgar_bmi_input_params->precipitation_mm_per_h * mm_to_m;
    bmi_unit_conv.volQ_timestep_m       = state->lgar_bmi_input_params->precipitation_mm_per_h * mm_to_m;
    bmi_unit_conv.volQ_gw_timestep_m    = 0.0;
    bmi_unit_conv.volPET_timestep_m     = 0.0;
    bmi_unit_conv.volrunoff_giuh_timestep_m = 0.0;
    bmi_unit_conv.volrunoff_giuh_ponded_m = 0.0;

    return;
  }

  // if lasam is coupled to soil freeze-thaw, frozen fraction module is called
  if (state->lgar_bmi_params.sft_coupled)
    frozen_factor_hydraulic_conductivity(state->lgar_bmi_params);

  double volchange_calib_cm = 0.0;

  if(state->lgar_bmi_params.calib_params_flag) {
    volchange_calib_cm = update_calibratable_parameters(); // change in soil water volume due to calibratable parameters
    state->lgar_bmi_params.calib_params_flag = false;
  }

  // local variables for readibility
  int subcycles;
  int num_layers = state->lgar_bmi_params.num_layers;

  // local variables for a full timestep (i.e., timestep of the forcing data)
  // see 'struct lgar_mass_balance_variables' in all.hxx for full description of the variables
  double precip_timestep_cm = 0.0;
  double PET_timestep_cm    = 0.0;
  double AET_timestep_cm    = 0.0;
  double volend_timestep_cm = lgar_calc_mass_bal(state->lgar_bmi_params.cum_layer_thickness_cm, state->head); // this should not be reset to 0.0 in the for loop
  double volin_timestep_cm  = 0.0;
  double volon_timestep_cm  = state->lgar_mass_balance.volon_timestep_cm;
  double volrunoff_timestep_cm      = 0.0;
  double volrech_timestep_cm        = 0.0;
  double surface_runoff_timestep_cm = 0.0; // direct surface runoff
  double volrunoff_giuh_timestep_cm = 0.0;
  double volrunoff_giuh_ponded_cm   = 0.0;
  double volQ_timestep_cm           = 0.0;
  double volQ_gw_timestep_cm        = 0.0;

  // local variables for a subtimestep (i.e., timestep of the model)
  double precip_subtimestep_cm;
  double precip_subtimestep_cm_per_h;
  double PET_subtimestep_cm;
  double PET_subtimestep_cm_per_h;
  double ponded_depth_subtimestep_cm;
  double AET_subtimestep_cm;
  double volstart_subtimestep_cm;
  double volend_subtimestep_cm = volend_timestep_cm; // this should not be reset to 0.0 in the for loop
  double volin_subtimestep_cm;
  double volon_subtimestep_cm;
  double volrunoff_subtimestep_cm;
  double volrech_subtimestep_cm;
  double surface_runoff_subtimestep_cm; // direct surface runoff
  double precip_previous_subtimestep_cm;
  double volQ_gw_subtimestep_cm = 0.0; // fix it for non-zero values after adding groundwater reservoir

  double subtimestep_h = state->lgar_bmi_params.timestep_h;
  int nint = state->lgar_bmi_params.nint;
  double wilting_point_psi_cm = state->lgar_bmi_params.wilting_point_psi_cm;
  double field_capacity_psi_cm = state->lgar_bmi_params.field_capacity_psi_cm;
  bool use_closed_form_G = state->lgar_bmi_params.use_closed_form_G;
  bool adaptive_timestep = state->lgar_bmi_params.adaptive_timestep;

  // constant value used in the AET function
  double AET_thresh_Theta = 0.85;    // scaled soil moisture (0-1) above which AET=PET (fix later!)
  double AET_expon        = 1.0;     // exponent that allows curvature of the rising portion of the Budyko curve (fix later!)
  double ponded_depth_max_cm = state->lgar_bmi_params.ponded_depth_max_cm;

  if (verbosity.compare("high") == 0) {
    bmilgar_ss <<"Pr  [cm/h] (timestep) = "<<state->lgar_bmi_input_params->precipitation_mm_per_h * mm_to_cm <<"\n";
    bmilgar_ss <<"PET [cm/h] (timestep) = "<<state->lgar_bmi_input_params->PET_mm_per_h * mm_to_cm <<"\n";
    LOG(bmilgar_ss.str(), LogLevel::INFO); bmilgar_ss.str("");
  }

  if (state->lgar_bmi_input_params->precipitation_mm_per_h < 0.0) {
    std::stringstream error_message;
    error_message << "Pr [mm/h] is less than 0: " << state->lgar_bmi_input_params->precipitation_mm_per_h;
    LOG(error_message.str(), LogLevel::INFO);
    throw std::runtime_error(error_message.str());
  }
  if (state->lgar_bmi_input_params->PET_mm_per_h < 0.0) {
    std::stringstream error_message;
    error_message << "PET [mm/h] is less than 0: " << state->lgar_bmi_input_params->PET_mm_per_h;
    LOG(error_message.str(), LogLevel::INFO);
    throw std::runtime_error(error_message.str());
  }

  // adaptive time step is set
  if (adaptive_timestep) {
    subtimestep_h = state->lgar_bmi_params.forcing_resolution_h;
    if (state->lgar_bmi_input_params->precipitation_mm_per_h > 10.0 || volon_timestep_cm > 0.0 ) {
      subtimestep_h = state->lgar_bmi_params.minimum_timestep_h;  //case where precip > 1 cm/h, or there is ponded head from the last time step
    }
    else if (state->lgar_bmi_input_params->precipitation_mm_per_h > 0.0) {
      subtimestep_h = state->lgar_bmi_params.minimum_timestep_h * 2.0;  //case where precip is less than 1 cm/h but greater than 0, and there is no ponded head
    }
    subtimestep_h = fmin(subtimestep_h, state->lgar_bmi_params.forcing_resolution_h);  //just in case the user has specified a minimum time step that would make the subtimestep_h greater than the forcing resolution
    state->lgar_bmi_params.timestep_h = subtimestep_h;
  }

  state->lgar_bmi_params.forcing_interval = int(state->lgar_bmi_params.forcing_resolution_h/state->lgar_bmi_params.timestep_h+1.0e-08); // add 1.0e-08 to prevent truncation error
  subcycles = state->lgar_bmi_params.forcing_interval;

  if (verbosity.compare("high") == 0) {
    LOG(LogLevel::DEBUG,"time step size in hours: %lf \n", state->lgar_bmi_params.timestep_h);
  }

  // subcycling loop (loop over model's timestep)
  for (int cycle=1; cycle <= subcycles; cycle++) {

    this->state->lgar_bmi_params.time_s    += subtimestep_h * state->units.hr_to_sec;
    this->state->lgar_bmi_params.timesteps ++;

    if (verbosity.compare("high") == 0 || verbosity.compare("low") == 0) {
      bmilgar_ss <<"BMI Update |---------------------------------------------------------------|\n";
      bmilgar_ss <<"BMI Update |Timesteps = "<< state->lgar_bmi_params.timesteps<<", Time [h] = "<<this->state->lgar_bmi_params.time_s / 3600.<<", Subcycle = "<< cycle <<" of "<<subcycles<<std::endl;
      LOG(bmilgar_ss.str(), LogLevel::INFO); bmilgar_ss.str("");
    }

    if( state->state_previous != NULL ){
      listDelete(state->state_previous);
      state->state_previous = NULL;
    }
    state->state_previous = listCopy(state->head);
    validate_wetting_front_state("after_state_copy", state);

    // ensure precip and PET are non-negative
    state->lgar_bmi_input_params->precipitation_mm_per_h = fmax(state->lgar_bmi_input_params->precipitation_mm_per_h, 0.0);
    state->lgar_bmi_input_params->PET_mm_per_h           = fmax(state->lgar_bmi_input_params->PET_mm_per_h, 0.0);

    precip_subtimestep_cm_per_h = state->lgar_bmi_input_params->precipitation_mm_per_h * mm_to_cm; // rate [cm/hour]
    PET_subtimestep_cm_per_h = state->lgar_bmi_input_params->PET_mm_per_h * mm_to_cm;

    ponded_depth_subtimestep_cm = precip_subtimestep_cm_per_h * subtimestep_h;
    ponded_depth_subtimestep_cm += volon_timestep_cm;

    precip_subtimestep_cm = precip_subtimestep_cm_per_h * subtimestep_h;
    PET_subtimestep_cm = PET_subtimestep_cm_per_h * subtimestep_h;

    if (verbosity.compare("high") == 0 || verbosity.compare("low") == 0) {
      bmilgar_ss <<"Pr [cm/h], Pr [cm] (subtimestep), subtimestep [h] = "<<state->lgar_bmi_input_params->precipitation_mm_per_h * mm_to_cm <<", "<< precip_subtimestep_cm <<", "<< subtimestep_h<<" ("<<subtimestep_h*3600<<" sec)"<<"\n";
      bmilgar_ss <<"PET [cm/h], PET [cm] (subtimestep) = "<<state->lgar_bmi_input_params->PET_mm_per_h * mm_to_cm <<", "<< PET_subtimestep_cm<<"\n";
      LOG(bmilgar_ss.str(), LogLevel::INFO); bmilgar_ss.str("");
    }

    AET_subtimestep_cm            = 0.0;
    volstart_subtimestep_cm       = 0.0;
    volin_subtimestep_cm          = 0.0;
    volrunoff_subtimestep_cm      = 0.0;
    volrech_subtimestep_cm        = 0.0;
    surface_runoff_subtimestep_cm = 0.0;

    precip_previous_subtimestep_cm = state->lgar_bmi_params.precip_previous_timestep_cm;

    num_layers = state->lgar_bmi_params.num_layers;
    double delta_theta;
    double dry_depth;

    if (PET_subtimestep_cm_per_h > 0.0) {
      AET_subtimestep_cm = calc_aet(PET_subtimestep_cm_per_h, subtimestep_h, wilting_point_psi_cm, field_capacity_psi_cm,
                                    state->lgar_bmi_params.layer_soil_type, AET_thresh_Theta, AET_expon,
                                    state->head, state->soil_properties);
    }

    precip_timestep_cm += precip_subtimestep_cm;
    PET_timestep_cm += fmax(PET_subtimestep_cm,0.0);

    volstart_subtimestep_cm = lgar_calc_mass_bal(state->lgar_bmi_params.cum_layer_thickness_cm, state->head);

    volon_timestep_cm = fmax(volon_timestep_cm,0.0);
    volon_timestep_cm = volon_timestep_cm > 1.0E-12 ? volon_timestep_cm : 0.0;

    int wf_free_drainage_demand = wetting_front_free_drainage(state->head);

    int soil_num = state->lgar_bmi_params.layer_soil_type[state->head->layer_num];
    double theta_e = state->soil_properties[soil_num].theta_e;
    bool is_top_wf_saturated = (state->head->theta+1.0E-12) >= theta_e ? true : false;

    bool create_surficial_front = (precip_previous_subtimestep_cm == 0.0 && precip_subtimestep_cm > 0.0);

    if (is_top_wf_saturated || volon_timestep_cm > 0.0)
      create_surficial_front = false;

    if (verbosity.compare("high") == 0 || verbosity.compare("low") == 0) {
      std::string flag        = (create_surficial_front && !is_top_wf_saturated) == true ? "Yes" : "No";
      std::string flag_top_wf = is_top_wf_saturated == true ? "Yes" : "No";
      bmilgar_ss <<"Is top wetting front saturated? "<< flag_top_wf  << "\n";
      bmilgar_ss <<"Create superficial wetting front? "<< flag << "\n";
      LOG(bmilgar_ss.str(), LogLevel::INFO); bmilgar_ss.str("");
    }

    if(create_surficial_front) {

      double temp_pd = 0.0;

      lgar_move_wetting_fronts(subtimestep_h, &temp_pd, wf_free_drainage_demand, volend_subtimestep_cm,
                               num_layers, &AET_subtimestep_cm, state->lgar_bmi_params.cum_layer_thickness_cm,
                               state->lgar_bmi_params.layer_soil_type, state->lgar_bmi_params.frozen_factor,
                               &state->head, state->state_previous, state->soil_properties);
      validate_wetting_front_state("after_move_wetting_fronts_create_path", state);

      if (temp_pd != 0.0){
        volrech_subtimestep_cm = temp_pd;
        volrech_timestep_cm += volrech_subtimestep_cm;
        temp_pd = 0.0;
      }

      dry_depth = lgar_calc_dry_depth(use_closed_form_G, nint, subtimestep_h, &delta_theta, state->lgar_bmi_params.layer_soil_type,
                                      state->lgar_bmi_params.cum_layer_thickness_cm, state->lgar_bmi_params.frozen_factor,
                                      state->head, state->soil_properties);

      if (verbosity.compare("high") == 0) {
        LOG(LogLevel::DEBUG,"State before moving creating new WF...\n");
        listPrint(state->head);
      }

      lgar_create_surficial_front(num_layers, &ponded_depth_subtimestep_cm, &volin_subtimestep_cm, dry_depth, state->head->theta,
                                  state->lgar_bmi_params.layer_soil_type, state->lgar_bmi_params.cum_layer_thickness_cm,
                                  state->lgar_bmi_params.frozen_factor, &state->head, state->soil_properties);
      validate_wetting_front_state("after_create_surficial_front", state);

      if (verbosity.compare("high") == 0) {
        LOG(LogLevel::DEBUG,"State after moving creating new WF...\n");
        listPrint(state->head);
      }

      if(state->state_previous != NULL ){
        listDelete(state->state_previous);
        state->state_previous = NULL;
      }
      state->state_previous = listCopy(state->head);
      validate_wetting_front_state("after_state_copy_post_create_surficial_front", state);

      volin_timestep_cm += volin_subtimestep_cm;

      if (verbosity.compare("high") == 0) {
        LOG(LogLevel::DEBUG,"New wetting front created...\n");
        listPrint(state->head);
      }
    }

    if (ponded_depth_subtimestep_cm > 0 && !create_surficial_front) {

      volrunoff_subtimestep_cm = lgar_insert_water(use_closed_form_G, nint, subtimestep_h, AET_subtimestep_cm, &ponded_depth_subtimestep_cm,
                                                   &volin_subtimestep_cm, precip_subtimestep_cm_per_h,
                                                   wf_free_drainage_demand, num_layers,
                                                   ponded_depth_max_cm, state->lgar_bmi_params.layer_soil_type,
                                                   state->lgar_bmi_params.cum_layer_thickness_cm,
                                                   state->lgar_bmi_params.frozen_factor, state->head,
                                                   state->soil_properties);
      validate_wetting_front_state("after_insert_water", state);

      volin_timestep_cm += volin_subtimestep_cm;
      volrunoff_timestep_cm += volrunoff_subtimestep_cm;
      volrech_subtimestep_cm = volin_subtimestep_cm;

      volon_subtimestep_cm = ponded_depth_subtimestep_cm;
      if (volrunoff_subtimestep_cm < 0) abort();
    }
    else {

      if (ponded_depth_subtimestep_cm < ponded_depth_max_cm) {
        volrunoff_timestep_cm += 0.0;
        volon_subtimestep_cm = ponded_depth_subtimestep_cm;
        ponded_depth_subtimestep_cm = 0.0;
        volrunoff_subtimestep_cm = 0.0;
      }
      else {
        volrunoff_subtimestep_cm = (ponded_depth_subtimestep_cm - ponded_depth_max_cm);
        volrunoff_timestep_cm += (ponded_depth_subtimestep_cm - ponded_depth_max_cm);
        volon_subtimestep_cm = ponded_depth_max_cm;
        ponded_depth_subtimestep_cm = ponded_depth_max_cm;
      }
    }

    if (!create_surficial_front) {
      double volin_subtimestep_cm_temp = volin_subtimestep_cm;
      lgar_move_wetting_fronts(subtimestep_h, &volin_subtimestep_cm, wf_free_drainage_demand, volend_subtimestep_cm,
                               num_layers, &AET_subtimestep_cm, state->lgar_bmi_params.cum_layer_thickness_cm,
                               state->lgar_bmi_params.layer_soil_type, state->lgar_bmi_params.frozen_factor,
                               &state->head, state->state_previous, state->soil_properties);
      validate_wetting_front_state("after_move_wetting_fronts", state);

      volrech_subtimestep_cm = volin_subtimestep_cm;
      volrech_timestep_cm += volrech_subtimestep_cm;

      volin_subtimestep_cm = volin_subtimestep_cm_temp;
    }

    lgar_dzdt_calc(use_closed_form_G, nint, ponded_depth_subtimestep_cm, state->lgar_bmi_params.layer_soil_type,
                   state->lgar_bmi_params.cum_layer_thickness_cm, state->lgar_bmi_params.frozen_factor,
                   state->head, state->soil_properties);
    validate_wetting_front_state("after_dzdt_calc", state);

    volend_subtimestep_cm = lgar_calc_mass_bal(state->lgar_bmi_params.cum_layer_thickness_cm, state->head);
    volend_timestep_cm = volend_subtimestep_cm;
    state->lgar_bmi_params.precip_previous_timestep_cm = precip_subtimestep_cm;

    double local_mb = volstart_subtimestep_cm + precip_subtimestep_cm + volon_timestep_cm - volrunoff_subtimestep_cm
                      - AET_subtimestep_cm - volon_subtimestep_cm - volrech_subtimestep_cm - volend_subtimestep_cm;

    AET_timestep_cm += AET_subtimestep_cm;
    volon_timestep_cm = volon_subtimestep_cm;

    surface_runoff_subtimestep_cm = volrunoff_subtimestep_cm;
    surface_runoff_timestep_cm += surface_runoff_subtimestep_cm ;

    volQ_gw_timestep_cm += volQ_gw_subtimestep_cm;

    if (verbosity.compare("high") == 0 || verbosity.compare("low") == 0) {
      LOG(LogLevel::DEBUG,"Printing wetting fronts at this subtimestep... \n");
      listPrint(state->head);
    }

    bool unexpected_local_error = fabs(local_mb) > 1.0E-4 ? true : false;

    if (verbosity.compare("high") == 0 || verbosity.compare("low") == 0 || unexpected_local_error) {
      LOG(LogLevel::DEBUG,"\nLocal mass balance at this timestep... \n\
      Error         = %14.10f \n\
      Initial water = %14.10f \n\
      Water added   = %14.10f \n\
      Ponded water  = %14.10f \n\
      Infiltration  = %14.10f \n\
      Runoff        = %14.10f \n\
      AET           = %14.10f \n\
      Percolation   = %14.10f \n\
      Final water   = %14.10f \n", local_mb, volstart_subtimestep_cm, precip_subtimestep_cm, volon_subtimestep_cm,
           volin_subtimestep_cm, volrunoff_subtimestep_cm, AET_subtimestep_cm, volrech_subtimestep_cm,
           volend_subtimestep_cm);

      if (unexpected_local_error) {
        LOG(LogLevel::DEBUG,"Local mass balance (in this timestep) is %14.10f, larger than expected, needs some debugging...\n ",local_mb);
        abort();
      }
    }

    state->lgar_mass_balance.local_mass_balance = local_mb;

    if (state->head->depth_cm <= 0.0) {
      std::stringstream error_message;
      error_message << "Cycle " << cycle << " has a depth less than or equal to 0: " << state->head->depth_cm;
      LOG(error_message.str(), LogLevel::INFO);
      throw std::runtime_error(error_message.str());
    }

    bool lasam_standalone = true;
#ifdef NGEN
    lasam_standalone = false;
#endif
    if ( (this->state->lgar_bmi_params.time_s >= this->state->lgar_bmi_params.endtime_s) && lasam_standalone)
      break;

  } // end of subcycling

  volrunoff_giuh_timestep_cm = giuh_convolution_integral(volrunoff_timestep_cm, num_giuh_ordinates, giuh_ordinates, giuh_runoff_queue);
  volQ_timestep_cm = volrunoff_giuh_timestep_cm;

  for (int i = 0; i < num_giuh_ordinates; ++i) {
    volrunoff_giuh_ponded_cm += giuh_runoff_queue[i];
  }

  validate_wetting_front_state("before_bmi_export", state);

  state->lgar_bmi_params.num_wetting_fronts = listLength(state->head);
  realloc_soil();

  struct wetting_front *current = state->head;
  for (int i=0; i<state->lgar_bmi_params.num_wetting_fronts; i++) {
    if (current == NULL) {
      std::stringstream error_message;
      error_message << "Wetting front at index " << i << " is null.";
      LOG(error_message.str(), LogLevel::INFO);
      throw std::runtime_error(error_message.str());
    }

    double theta_export = current->theta;
    double depth_m_export = current->depth_cm * state->units.cm_to_m;

    if (!FIRST_BAD_FOUND &&
        (!std::isfinite(theta_export) || theta_export > 1.0 || theta_export < 0.0)) {

      FIRST_BAD_FOUND = true;

      std::stringstream msg;
      msg << "FIRST FAILURE DETECTED"
    << " timestep=" << state->lgar_bmi_params.timesteps
    << " wf_index=" << i
    << " theta=" << theta_export
    << " psi_cm=" << current->psi_cm
    << " depth_cm=" << current->depth_cm
    << " layer=" << current->layer_num;

      LOG(msg.str(), LogLevel::SEVERE);
      abort();
    }

    static int bad_export_log_count = 0;
    static const int bad_export_log_limit = 20;

    if ((!std::isfinite(theta_export) || theta_export < 0.0 || theta_export > 1.0) &&
        bad_export_log_count < bad_export_log_limit) {
      bad_export_log_count++;

      std::stringstream msg;
      msg << "BAD EXPORT"
          << " wf=" << i
          << " theta=" << theta_export
          << " psi_cm=" << current->psi_cm
          << " depth_cm=" << current->depth_cm
          << " layer=" << current->layer_num
          << " dzdt=" << current->dzdt_cm_per_h
          << " to_bottom=" << current->to_bottom
          << " timesteps=" << state->lgar_bmi_params.timesteps;
      LOG(msg.str(), LogLevel::INFO);
    }

    state->lgar_bmi_params.soil_moisture_wetting_fronts[i] = theta_export;
    state->lgar_bmi_params.soil_depth_wetting_fronts[i] = depth_m_export;

    if (!std::isfinite(theta_export) || theta_export > 1.0 || theta_export < 0.0) {
      if (GLOBAL_ERROR_COUNT < GLOBAL_ERROR_LIMIT) {
        GLOBAL_ERROR_COUNT++;

        std::stringstream msg;
        msg << "CRITICAL: INVALID THETA DETECTED"
            << " timestep=" << state->lgar_bmi_params.timesteps
            << " wf_index=" << i
            << " theta=" << theta_export
            << " psi_cm=" << current->psi_cm
            << " depth_cm=" << current->depth_cm
            << " layer=" << current->layer_num;
        LOG(msg.str(), LogLevel::SEVERE);
      }
    }

    if (!std::isfinite(state->lgar_bmi_params.soil_moisture_wetting_fronts[i])) {
      std::stringstream msg;
      msg << "LASAM export buffer became non-finite immediately after write"
          << " i=" << i
          << " theta_export=" << theta_export
          << " timestep=" << state->lgar_bmi_params.timesteps;
      LOG(msg.str(), LogLevel::INFO);
    }

    int layer_num = current->layer_num;
    int soil_num_local = -1;
    double theta_e_local = -9999.0;
    double theta_r_local = -9999.0;

    if (layer_num >= 1 && layer_num <= state->lgar_bmi_params.num_layers) {
      soil_num_local = state->lgar_bmi_params.layer_soil_type[layer_num];
      if (soil_num_local >= 1 && soil_num_local <= state->lgar_bmi_params.num_soil_types) {
        theta_e_local = state->soil_properties[soil_num_local].theta_e;
        theta_r_local = state->soil_properties[soil_num_local].theta_r;
      }
    }

    bool suspicious = false;

    if (!std::isfinite(theta_export) || !std::isfinite(depth_m_export) || !std::isfinite(current->psi_cm)) {
      suspicious = true;
    }
    if (theta_export < 0.0) {
      suspicious = true;
    }
    if (theta_e_local > 0.0 && theta_export > theta_e_local + 1.0e-10) {
      suspicious = true;
    }
    if (current->depth_cm < 0.0) {
      suspicious = true;
    }

    if (suspicious) {
      std::stringstream msg;
      msg << "LASAM suspicious wf export:"
          << " timestep=" << state->lgar_bmi_params.timesteps
          << " wf_index=" << i
          << " layer=" << current->layer_num
          << " soil=" << soil_num_local
          << " theta=" << theta_export
          << " theta_r=" << theta_r_local
          << " theta_e=" << theta_e_local
          << " psi_cm=" << current->psi_cm
          << " depth_cm=" << current->depth_cm
          << " depth_m=" << depth_m_export
          << " dzdt_cm_per_h=" << current->dzdt_cm_per_h
          << " to_bottom=" << current->to_bottom;
      LOG(msg.str(), LogLevel::INFO);

      struct wetting_front *prev = state->state_previous;
      struct wetting_front *prev_match = NULL;
      int prev_idx = 0;
      while (prev != NULL) {
        if (prev_idx == i) {
          prev_match = prev;
          break;
        }
        prev = prev->next;
        prev_idx++;
      }

      if (prev_match != NULL) {
        std::stringstream msg_prev;
        msg_prev << "LASAM previous wf state:"
                 << " timestep=" << (state->lgar_bmi_params.timesteps - 1)
                 << " wf_index=" << i
                 << " layer=" << prev_match->layer_num
                 << " theta=" << prev_match->theta
                 << " psi_cm=" << prev_match->psi_cm
                 << " depth_cm=" << prev_match->depth_cm
                 << " dzdt_cm_per_h=" << prev_match->dzdt_cm_per_h
                 << " to_bottom=" << prev_match->to_bottom;
        LOG(msg_prev.str(), LogLevel::INFO);
      }
    }

    if (verbosity.compare("high") == 0) {
      bmilgar_ss <<"Wetting fronts (bmi outputs) (depth in meters, theta)= "
                 <<state->lgar_bmi_params.soil_depth_wetting_fronts[i]
                 <<" "<<state->lgar_bmi_params.soil_moisture_wetting_fronts[i]<<"\n";
      LOG(bmilgar_ss.str(), LogLevel::INFO); bmilgar_ss.str("");
    }

    current = current->next;
  }

  state->lgar_mass_balance.volprecip_timestep_cm  = precip_timestep_cm;
  state->lgar_mass_balance.volin_timestep_cm      = volin_timestep_cm;
  state->lgar_mass_balance.volon_timestep_cm      = volon_timestep_cm;
  state->lgar_mass_balance.volend_timestep_cm     = volend_timestep_cm;
  state->lgar_mass_balance.volAET_timestep_cm     = AET_timestep_cm;
  state->lgar_mass_balance.volrech_timestep_cm    = volrech_timestep_cm;
  state->lgar_mass_balance.volrunoff_timestep_cm  = volrunoff_timestep_cm;
  state->lgar_mass_balance.volQ_timestep_cm       = volQ_timestep_cm;
  state->lgar_mass_balance.volQ_gw_timestep_cm    = volQ_gw_timestep_cm;
  state->lgar_mass_balance.volPET_timestep_cm     = PET_timestep_cm;
  state->lgar_mass_balance.volrunoff_giuh_timestep_cm = volrunoff_giuh_timestep_cm;

  state->lgar_mass_balance.volprecip_cm  += precip_timestep_cm;
  state->lgar_mass_balance.volin_cm      += volin_timestep_cm;
  state->lgar_mass_balance.volon_cm       = volon_timestep_cm;
  state->lgar_mass_balance.volend_cm      = volend_timestep_cm;
  state->lgar_mass_balance.volAET_cm     += AET_timestep_cm;
  state->lgar_mass_balance.volrech_cm    += volrech_timestep_cm;
  state->lgar_mass_balance.volrunoff_cm  += volrunoff_timestep_cm;
  state->lgar_mass_balance.volQ_cm       += volQ_timestep_cm;
  state->lgar_mass_balance.volQ_gw_cm    += volQ_gw_timestep_cm;
  state->lgar_mass_balance.volPET_cm     += PET_timestep_cm;
  state->lgar_mass_balance.volrunoff_giuh_cm  += volrunoff_giuh_timestep_cm;
  state->lgar_mass_balance.volchange_calib_cm += volchange_calib_cm ;

  bmi_unit_conv.mass_balance_m        = state->lgar_mass_balance.local_mass_balance * state->units.cm_to_m;
  bmi_unit_conv.volprecip_timestep_m  = precip_timestep_cm * state->units.cm_to_m;
  bmi_unit_conv.volin_timestep_m      = volin_timestep_cm * state->units.cm_to_m;
  bmi_unit_conv.volend_timestep_m     = volend_timestep_cm * state->units.cm_to_m;
  bmi_unit_conv.volAET_timestep_m     = AET_timestep_cm * state->units.cm_to_m;
  bmi_unit_conv.volrech_timestep_m    = volrech_timestep_cm * state->units.cm_to_m;
  bmi_unit_conv.volrunoff_timestep_m  = volrunoff_timestep_cm * state->units.cm_to_m;
  bmi_unit_conv.volQ_timestep_m       = volQ_timestep_cm * state->units.cm_to_m;
  bmi_unit_conv.volQ_gw_timestep_m    = volQ_gw_timestep_cm * state->units.cm_to_m;
  bmi_unit_conv.volPET_timestep_m     = PET_timestep_cm * state->units.cm_to_m;
  bmi_unit_conv.volrunoff_giuh_timestep_m = volrunoff_giuh_timestep_cm * state->units.cm_to_m;
  bmi_unit_conv.volrunoff_giuh_ponded_m = volrunoff_giuh_ponded_cm * state->units.cm_to_m;
}

void BmiLGAR::
UpdateUntil(double t)
{
  if (t <= 0.0) {
    const char *error_message = "Time must be greater than 0.";
    LOG(LogLevel::INFO, error_message);
    throw std::invalid_argument(error_message);
  }
  this->Update();
}

struct model_state* BmiLGAR::get_model()
{
  return state;
}

void BmiLGAR::
global_mass_balance()
{
  lgar_global_mass_balance(this->state, giuh_runoff_queue);
}

double BmiLGAR::
update_calibratable_parameters()
{
  int soil, layer_num;
  struct wetting_front *current = state->head;

  if (verbosity.compare("high") == 0)
    listPrint(state->head);
  
  double volstart_before = lgar_calc_mass_bal(state->lgar_bmi_params.cum_layer_thickness_cm, state->head);

  for (int i=0; i<state->lgar_bmi_params.num_wetting_fronts; i++) {//first we update the parameters that depend on soil layer, for each layer
    layer_num  = current->layer_num;
    soil = state->lgar_bmi_params.layer_soil_type[layer_num];
    
    if (current == NULL) {
      std::stringstream error_message;
      error_message << "Wetting front at index " << i << " is null.";
      LOG(error_message.str(), LogLevel::INFO);
      throw std::invalid_argument(error_message.str());
    }

    if (verbosity.compare("high") == 0 || verbosity.compare("low") == 0) {
      bmilgar_ss <<"----------- Calibratable parameters depending on soil layer (initial values) ----------- \n";
      bmilgar_ss <<"| soil_type = "<< soil <<", layer = "<<layer_num
	       <<", smcmax = "   << state->soil_properties[soil].theta_e
	       <<", smcmin = "   << state->soil_properties[soil].theta_r
	       <<", vg_n = "     << state->soil_properties[soil].vg_n
	       <<", vg_alpha = " << state->soil_properties[soil].vg_alpha_per_cm
	       <<", Ksat = "     << state->soil_properties[soil].Ksat_cm_per_h
	       <<", theta = "    << current->theta <<"\n";
        LOG(bmilgar_ss.str(), LogLevel::INFO); bmilgar_ss.str("");  
    }
    
    state->soil_properties[soil].theta_e = state->lgar_calib_params.theta_e[layer_num-1];
    state->soil_properties[soil].theta_r = state->lgar_calib_params.theta_r[layer_num-1];
    state->soil_properties[soil].vg_n    = state->lgar_calib_params.vg_n[layer_num-1];
    state->soil_properties[soil].vg_m    = 1.0 - 1.0/state->soil_properties[soil].vg_n;
    state->soil_properties[soil].vg_alpha_per_cm = state->lgar_calib_params.vg_alpha[layer_num-1];
    state->soil_properties[soil].Ksat_cm_per_h   = state->lgar_calib_params.Ksat[layer_num-1];
    
    current->theta = calc_theta_from_h(current->psi_cm, state->soil_properties[soil].vg_alpha_per_cm,
				       state->soil_properties[soil].vg_m, state->soil_properties[soil].vg_n,
				       state->soil_properties[soil].theta_e, state->soil_properties[soil].theta_r);

    if (verbosity.compare("high") == 0 || verbosity.compare("low") == 0) {
      bmilgar_ss <<"----------- Calibratable parameters depending on soil layer (updated values) ----------- \n";
      bmilgar_ss <<"| soil_type = "<< soil <<", layer = "<<layer_num
	       <<", smcmax = "   << state->soil_properties[soil].theta_e
	       <<", smcmin = "   << state->soil_properties[soil].theta_r
	       <<", vg_n = "     << state->soil_properties[soil].vg_n
	       <<", vg_alpha = " << state->soil_properties[soil].vg_alpha_per_cm
	       <<", Ksat = "     << state->soil_properties[soil].Ksat_cm_per_h
	       <<", theta = "    << current->theta <<"\n";
        LOG(bmilgar_ss.str(), LogLevel::INFO); bmilgar_ss.str("");  
    }
    
    current = current->next;
  }

  //next we update the parameters that apply to the whole model domain and do not depend on soil layer
  if (verbosity.compare("high") == 0 || verbosity.compare("low") == 0) {
    bmilgar_ss <<"----------- Calibratable parameters independent of soil layer (initial values) ----------- \n";
    bmilgar_ss <<"field_capacity_psi = "   << state->lgar_bmi_params.field_capacity_psi_cm
      <<", ponded_depth_max = "     << state->lgar_bmi_params.ponded_depth_max_cm <<"\n";
    LOG(bmilgar_ss.str(), LogLevel::INFO); bmilgar_ss.str("");  
  }

  state->lgar_bmi_params.field_capacity_psi_cm = state->lgar_calib_params.field_capacity_psi;
  state->lgar_bmi_params.ponded_depth_max_cm   = state->lgar_calib_params.ponded_depth_max;

  if (verbosity.compare("high") == 0 || verbosity.compare("low") == 0) {
    bmilgar_ss <<"----------- Calibratable parameters independent of soil layer (updated values) ----------- \n";
    bmilgar_ss <<"field_capacity_psi = "   << state->lgar_bmi_params.field_capacity_psi_cm
      <<", ponded_depth_max = "     << state->lgar_bmi_params.ponded_depth_max_cm <<"\n";
    LOG(bmilgar_ss.str(), LogLevel::INFO); bmilgar_ss.str("");  
  }
  
  if (verbosity.compare("high") == 0)
    listPrint(state->head);
  
  double volstart_after = lgar_calc_mass_bal(state->lgar_bmi_params.cum_layer_thickness_cm, state->head);

  if (verbosity.compare("high") == 0 || verbosity.compare("low") == 0) {
    bmilgar_ss <<"Mass of water (before and after) = "<< volstart_before<<", "<< volstart_after <<"\n";
    LOG(bmilgar_ss.str(), LogLevel::INFO); bmilgar_ss.str("");  
  }
  
  return volstart_after - volstart_before;
}

void BmiLGAR::
Finalize()
{
  global_mass_balance();
  listDelete(state->head);
  listDelete(state->state_previous);

  delete [] state->soil_properties;

  if (state->lgar_bmi_params.soil_depth_wetting_fronts != nullptr) {
    delete [] state->lgar_bmi_params.soil_depth_wetting_fronts;
    state->lgar_bmi_params.soil_depth_wetting_fronts = nullptr;
  }

  if (state->lgar_bmi_params.soil_moisture_wetting_fronts != nullptr) {
    delete [] state->lgar_bmi_params.soil_moisture_wetting_fronts;
    state->lgar_bmi_params.soil_moisture_wetting_fronts = nullptr;
  }

  g_lasam_export_capacity.erase(state);

  delete [] state->lgar_bmi_params.soil_temperature;
  delete [] state->lgar_bmi_params.soil_temperature_z;
  delete [] state->lgar_bmi_params.layer_soil_type;

  delete [] state->lgar_calib_params.theta_e;
  delete [] state->lgar_calib_params.theta_r;
  delete [] state->lgar_calib_params.vg_n;
  delete [] state->lgar_calib_params.vg_alpha;
  delete [] state->lgar_calib_params.Ksat;

  delete [] state->lgar_bmi_params.layer_thickness_cm;
  delete [] state->lgar_bmi_params.cum_layer_thickness_cm;
  delete [] state->lgar_bmi_params.giuh_ordinates;
  delete [] state->lgar_bmi_params.frozen_factor;
  delete state->lgar_bmi_input_params;
  delete state;
  this->state = NULL;
}

int BmiLGAR::
GetVarGrid(std::string name)
{
  if (
    name.compare("soil_storage_model") == 0
    || name.compare("soil_num_wetting_fronts") == 0
    || name.compare("num_wetting_fronts") == 0
    || name.compare("serialization_free") == 0
  ) // int
    return 0;
  else if (
    name.compare("precipitation_rate") == 0
    || name.compare("precipitation_rate_out") == 0
    || name.compare("precipitation") == 0
    || name.compare("potential_evapotranspiration_rate") == 0
    || name.compare("potential_evapotranspiration") == 0
    || name.compare("actual_evapotranspiration") == 0
    || name.compare("surface_runoff") == 0
    || name.compare("giuh_runoff") == 0
    || name.compare("soil_storage") == 0
    || name.compare("field_capacity") == 0
    || name.compare("ponded_depth_max") == 0
    || name.compare("total_discharge") == 0
    || name.compare("infiltration") == 0
    || name.compare("percolation") == 0
    || name.compare("groundwater_to_stream_recharge") == 0
    || name.compare("mass_balance") == 0
    || name.compare(NWM_PONDED_DEPTH_OUT_VAR) == 0
    || name.compare("reset_time") == 0
  ) // double
    return 1;
  else if (
    name.compare("soil_depth_layers") == 0
    || name.compare("smcmax") == 0
    || name.compare("smcmin") == 0
    || name.compare("van_genuchten_m") == 0
    || name.compare("van_genuchten_alpha") == 0
    || name.compare("van_genuchten_n") == 0
    || name.compare("hydraulic_conductivity") == 0
  ) // array of doubles (fixed length)
    return 2;
  else if (
    name.compare("soil_moisture_wetting_fronts") == 0
    || name.compare("soil_depth_wetting_fronts") == 0
  ) // array of doubles (dynamic length)
    return 3;
  else if (
    name.compare("soil_temperature_profile") == 0
  ) // array of doubles (fixed and of the size of soil temperature profile)
    return 4;
  else if (
    name.compare("serialization_state") == 0
  ) // char
    return 5;
  else if (
    name.compare("serialization_create") == 0
    || name.compare("serialization_size") == 0
  ) // uint64_t
    return 6;
  else
    return -1;
}


std::string BmiLGAR::
GetVarType(std::string name)
{
  int var_grid = GetVarGrid(name);

  if (var_grid == 0)
    return "int";
  else if (var_grid == 1 || var_grid == 2 || var_grid == 3 || var_grid == 4)
    return "double";
  else if (var_grid == 5)
    return "char";
  else if (var_grid == 6)
    return "uint64_t";
  else
    return "none";
}


int BmiLGAR::
GetVarItemsize(std::string name)
{
  int var_grid = GetVarGrid(name);

   if (var_grid == 0)
    return sizeof(int);
  else if (var_grid == 1 || var_grid == 2 || var_grid == 3 || var_grid == 4)
    return sizeof(double);
  else if (var_grid == 5)
    return sizeof(char);
  else if (var_grid == 6)
    return sizeof(uint64_t);
  else
    return 0;
}


std::string BmiLGAR::
GetVarUnits(std::string name)
{
  if (name.compare("precipitation_rate") == 0 || name.compare("precipitation_rate_out") == 0 
           || name.compare("potential_evapotranspiration_rate") == 0)
    return "mm h^-1";
  else if (name.compare("precipitation") == 0 || name.compare("potential_evapotranspiration") == 0
	   || name.compare("actual_evapotranspiration") == 0) // double
    return "m";
  else if (name.compare("surface_runoff") == 0 || name.compare("giuh_runoff") == 0
	   || name.compare("soil_storage") == 0 || name.compare(NWM_PONDED_DEPTH_OUT_VAR) == 0) // double
    return "m";
  else if (name.compare("total_discharge") == 0 || name.compare("infiltration") == 0
	   || name.compare("percolation") == 0) // double
    return "m";
  else if (name.compare("mass_balance") == 0 || name.compare("groundwater_to_stream_recharge") == 0)
    return "m";
  else if (name.compare("soil_moisture_wetting_fronts") == 0) // array of doubles
    return "none";
  else if (name.compare("soil_depth_layers") == 0 || name.compare("soil_depth_wetting_fronts") == 0) // array of doubles
    return "m";
  else if (name.compare("soil_temperature_profile") == 0)
    return "K";
  else
    return "none";

}


int BmiLGAR::
GetVarNbytes(std::string name)
{
  int itemsize;
  int gridsize;

  itemsize = this->GetVarItemsize(name);
  gridsize = this->GetGridSize(this->GetVarGrid(name));
  return itemsize * gridsize;
}


std::string BmiLGAR::
GetVarLocation(std::string name)
{
  if (name.compare("precipitation_rate") == 0 || name.compare("precipitation_rate_out") == 0 ||
    name.compare("precipitation") == 0 ||
    name.compare("potential_evapotranspiration") == 0 ||
    name.compare("potential_evapotranspiration_rate") == 0 ||
    name.compare("actual_evapotranspiration") == 0)
    return "node";
  else if (name.compare("surface_runoff") == 0 || name.compare("giuh_runoff") == 0
	   || name.compare("soil_storage") == 0 || name.compare(NWM_PONDED_DEPTH_OUT_VAR) == 0) // double
    return "node";
   else if (name.compare("total_discharge") == 0 || name.compare("infiltration") == 0
	    || name.compare("percolation") == 0 || name.compare("groundwater_to_stream_recharge") == 0) // double
    return "node";
  else if (name.compare("soil_moisture_wetting_fronts") == 0) // array of doubles
    return "node";
  else if (name.compare("mass_balance") == 0)
    return "node";
  else if (name.compare("soil_depth_layers") == 0 || name.compare("soil_depth_wetting_fronts") == 0
	   || name.compare("soil_num_wetting_fronts") == 0) // array of doubles
    return "node";
  else if (name.compare("soil_temperature_profile") == 0)
    return "node";
  else
    return "none";
}

void BmiLGAR::
GetGridShape(const int grid, int *shape)
{
  if (grid == 2)
    shape[0] = this->state->lgar_bmi_params.num_layers;
  else if (grid == 3) // number of wetting fronts (dynamic)
    shape[1] = this->state->lgar_bmi_params.num_wetting_fronts;
}


void BmiLGAR::
GetGridSpacing (const int grid, double * spacing)
{
  if (grid == 0) {
    spacing[0] = this->state->lgar_bmi_params.spacing[0];
  }
}


void BmiLGAR::
GetGridOrigin (const int grid, double *origin)
{
  if (grid == 0) {
    origin[0] = this->state->lgar_bmi_params.origin[0];
  }
}


int BmiLGAR::
GetGridRank(const int grid)
{
  if (grid == 0 || grid == 1 || grid == 2 || grid == 3 || grid == 4)
    return 1;
  else
    return -1;
}


int BmiLGAR::
GetGridSize(const int grid)
{
  if (grid == 0 || grid == 1 || grid == 6)
    return 1;
  else if (grid == 2) // number of layers (fixed)
    return this->state->lgar_bmi_params.num_layers;
  else if (grid == 3) // number of wetting fronts (dynamic)
    return this->state->lgar_bmi_params.num_wetting_fronts;
  else if (grid == 4) // number of cells (discretized temperature profile, input from SFT)
    return this->state->lgar_bmi_params.num_cells_temp;
  else if (grid == 5)
    return this->m_serialized_length; // serialized state length
  else
    return -1;
}

void BmiLGAR::
GetValue (std::string name, void *dest)
{
  if (dest == NULL) {
    std::stringstream errMsg;
    errMsg << "GetValue: destination pointer is null for variable " << name;
    LOG(errMsg.str(), LogLevel::INFO);
    throw std::runtime_error(errMsg.str());
  }

  // For dynamic wetting-front arrays, build directly from the linked-list state.
  // This avoids relying on exported buffers that may be stale/corrupted.
  if (name.compare("soil_moisture_wetting_fronts") == 0) {
    int count = state->lgar_bmi_params.num_wetting_fronts;
    double *out = static_cast<double*>(dest);
    struct wetting_front *current = state->head;

    for (int i = 0; i < count; i++) {
      if (current == NULL) {
        std::stringstream errMsg;
        errMsg << "GetValue: wetting front list ended early for variable " << name
               << " at index " << i << " of " << count;
        LOG(errMsg.str(), LogLevel::INFO);
        throw std::runtime_error(errMsg.str());
      }
      out[i] = current->theta;
      current = current->next;
    }
    return;
  }

  if (name.compare("soil_depth_wetting_fronts") == 0) {
    int count = state->lgar_bmi_params.num_wetting_fronts;
    double *out = static_cast<double*>(dest);
    struct wetting_front *current = state->head;

    for (int i = 0; i < count; i++) {
      if (current == NULL) {
        std::stringstream errMsg;
        errMsg << "GetValue: wetting front list ended early for variable " << name
               << " at index " << i << " of " << count;
        LOG(errMsg.str(), LogLevel::INFO);
        throw std::runtime_error(errMsg.str());
      }
      out[i] = current->depth_cm * state->units.cm_to_m;
      current = current->next;
    }
    return;
  }

  void * src = NULL;
  int nbytes = 0;

  src = this->GetValuePtr(name);
  nbytes = this->GetVarNbytes(name);

  if (src == NULL) {
    std::stringstream errMsg;
    errMsg << "GetValue: source pointer is null for variable " << name;
    LOG(errMsg.str(), LogLevel::INFO);
    throw std::runtime_error(errMsg.str());
  }

  memcpy(dest, src, nbytes);
}

void *BmiLGAR::
GetValuePtr (std::string name)
{
  if (name.compare("precipitation_rate") == 0)
    return (void*)(&this->state->lgar_bmi_input_params->precipitation_mm_per_h);
  else if (name.compare("precipitation_rate_out") == 0)
    return (void*)(&this->state->lgar_bmi_input_params->precipitation_mm_per_h);
  else if (name.compare("precipitation") == 0)
    return (void*)(&bmi_unit_conv.volprecip_timestep_m);
  else if (name.compare("potential_evapotranspiration_rate") == 0)
    return (void*)(&this->state->lgar_bmi_input_params->PET_mm_per_h);
  else if (name.compare("potential_evapotranspiration") == 0)
    return (void*)(&bmi_unit_conv.volPET_timestep_m);
  else if (name.compare("actual_evapotranspiration") == 0)
    return (void*)(&bmi_unit_conv.volAET_timestep_m);
  else if (name.compare("surface_runoff") == 0)
    return (void*)(&bmi_unit_conv.volrunoff_timestep_m);
  else if (name.compare("giuh_runoff") == 0)
    return (void*)(&bmi_unit_conv.volrunoff_giuh_timestep_m);
  else if (name.compare("soil_storage") == 0)
    return (void*)(&bmi_unit_conv.volend_timestep_m);
  else if (name.compare("total_discharge") == 0)
    return (void*)(&bmi_unit_conv.volQ_timestep_m);
  else if (name.compare("infiltration") == 0)
    return (void*)(&bmi_unit_conv.volin_timestep_m);
  else if (name.compare("percolation") == 0)
    return (void*)(&bmi_unit_conv.volrech_timestep_m);
  else if (name.compare("groundwater_to_stream_recharge") == 0)
    return (void*)(&bmi_unit_conv.volQ_gw_timestep_m);
  else if (name.compare("mass_balance") == 0)
    return (void*)(&bmi_unit_conv.mass_balance_m);
  else if (name.compare(NWM_PONDED_DEPTH_OUT_VAR) == 0)
    return (void*)(&this->bmi_unit_conv.volrunoff_giuh_ponded_m);
  else if (name.compare("soil_depth_layers") == 0)
    return (void*)this->state->lgar_bmi_params.cum_layer_thickness_cm;
  else if (name.compare("soil_moisture_wetting_fronts") == 0) {
    // Refresh export buffer from linked-list state before returning pointer.
    int count = state->lgar_bmi_params.num_wetting_fronts;
    this->realloc_soil();

    struct wetting_front *current = state->head;
    for (int i = 0; i < count; i++) {
      if (current == NULL) {
        std::stringstream errMsg;
        errMsg << "GetValuePtr: wetting front list ended early for variable "
               << name << " at index " << i << " of " << count;
        LOG(errMsg.str(), LogLevel::INFO);
        throw std::runtime_error(errMsg.str());
      }
      this->state->lgar_bmi_params.soil_moisture_wetting_fronts[i] = current->theta;
      current = current->next;
    }
    return (void*)this->state->lgar_bmi_params.soil_moisture_wetting_fronts;
  }
  else if (name.compare("soil_depth_wetting_fronts") == 0) {
    // Refresh export buffer from linked-list state before returning pointer.
    int count = state->lgar_bmi_params.num_wetting_fronts;
    this->realloc_soil();

    struct wetting_front *current = state->head;
    for (int i = 0; i < count; i++) {
      if (current == NULL) {
        std::stringstream errMsg;
        errMsg << "GetValuePtr: wetting front list ended early for variable "
               << name << " at index " << i << " of " << count;
        LOG(errMsg.str(), LogLevel::INFO);
        throw std::runtime_error(errMsg.str());
      }
      this->state->lgar_bmi_params.soil_depth_wetting_fronts[i] = current->depth_cm * state->units.cm_to_m;
      current = current->next;
    }
    return (void*)this->state->lgar_bmi_params.soil_depth_wetting_fronts;
  }
  else if (name.compare("soil_num_wetting_fronts") == 0)
    return (void*)(&state->lgar_bmi_params.num_wetting_fronts);
  else if (name.compare("num_wetting_fronts") == 0)
    return (void*)(&state->lgar_bmi_params.num_wetting_fronts);
  else if (name.compare("soil_temperature_profile") == 0)
    return (void*)this->state->lgar_bmi_params.soil_temperature;
  else if (name.compare("smcmax") == 0)
    return (void*)this->state->lgar_calib_params.theta_e;
  else if (name.compare("smcmin") == 0)
    return (void*)this->state->lgar_calib_params.theta_r;
  else if (name.compare("van_genuchten_n") == 0)
    return (void*)this->state->lgar_calib_params.vg_n;
  else if (name.compare("van_genuchten_alpha") == 0)
    return (void*)this->state->lgar_calib_params.vg_alpha;
  else if (name.compare("hydraulic_conductivity") == 0)
    return (void*)this->state->lgar_calib_params.Ksat;
  else if (name.compare("ponded_depth_max") == 0)
    return (void*)&this->state->lgar_calib_params.ponded_depth_max;
  else if (name.compare("field_capacity") == 0)
    return (void*)&this->state->lgar_calib_params.field_capacity_psi;
  else if (name.compare("serialization_state") == 0)
    return (void*)(this->m_serialized.data());
  else if (name.compare("serialization_size") == 0) {
    return (void*)(&this->m_serialized_length);
  } else {
    std::stringstream errMsg;
    errMsg << "variable "<< name << " does not exist";
    throw std::runtime_error(errMsg.str());
    return NULL;
  }

  return NULL;
}

void BmiLGAR::
GetValueAtIndices (std::string name, void *dest, int *inds, int len)
{
  if (dest == NULL) {
    std::stringstream errMsg;
    errMsg << "GetValueAtIndices: destination pointer is null for variable " << name;
    LOG(errMsg.str(), LogLevel::INFO);
    throw std::runtime_error(errMsg.str());
  }

  if (inds == NULL) {
    std::stringstream errMsg;
    errMsg << "GetValueAtIndices: indices pointer is null for variable " << name;
    LOG(errMsg.str(), LogLevel::INFO);
    throw std::runtime_error(errMsg.str());
  }

  if (len < 0) {
    std::stringstream errMsg;
    errMsg << "GetValueAtIndices: invalid len=" << len << " for variable " << name;
    LOG(errMsg.str(), LogLevel::INFO);
    throw std::runtime_error(errMsg.str());
  }

  int itemsize = this->GetVarItemsize(name);
  if (itemsize <= 0) {
    std::stringstream errMsg;
    errMsg << "GetValueAtIndices: invalid itemsize=" << itemsize
           << " for variable " << name;
    LOG(errMsg.str(), LogLevel::INFO);
    throw std::runtime_error(errMsg.str());
  }

  // For dynamic wetting-front arrays, rebuild values through GetValue() so we do not
  // rely on any potentially stale export-buffer pointer.
  if (name.compare("soil_moisture_wetting_fronts") == 0 ||
      name.compare("soil_depth_wetting_fronts") == 0) {

    int count = state->lgar_bmi_params.num_wetting_fronts;
    if (count <= 0) {
      std::stringstream errMsg;
      errMsg << "GetValueAtIndices: invalid num_wetting_fronts=" << count
             << " for variable " << name;
      LOG(errMsg.str(), LogLevel::INFO);
      throw std::runtime_error(errMsg.str());
    }

    std::vector<double> values(count);
    this->GetValue(name, values.data());

    char *ptr = (char *)dest;
    for (int i = 0; i < len; i++, ptr += itemsize) {
      if (inds[i] < 0 || inds[i] >= count) {
        std::stringstream errMsg;
        errMsg << "GetValueAtIndices: index " << inds[i]
               << " out of bounds [0," << (count - 1) << "] for variable " << name;
        LOG(errMsg.str(), LogLevel::INFO);
        throw std::runtime_error(errMsg.str());
      }

      memcpy(ptr, &values[inds[i]], itemsize);
    }
    return;
  }

  void * src = this->GetValuePtr(name);

  if (src == NULL) {
    std::stringstream errMsg;
    errMsg << "GetValueAtIndices: source pointer is null for variable " << name;
    LOG(errMsg.str(), LogLevel::INFO);
    throw std::runtime_error(errMsg.str());
  }

  int gridsize = this->GetGridSize(this->GetVarGrid(name));
  char *ptr = (char *)dest;

  for (int i = 0; i < len; i++, ptr += itemsize) {
    if (inds[i] < 0 || inds[i] >= gridsize) {
      std::stringstream errMsg;
      errMsg << "GetValueAtIndices: index " << inds[i]
             << " out of bounds [0," << (gridsize - 1) << "] for variable " << name;
      LOG(errMsg.str(), LogLevel::INFO);
      throw std::runtime_error(errMsg.str());
    }

    int offset = inds[i] * itemsize;
    memcpy(ptr, (char *)src + offset, itemsize);
  }
}

void BmiLGAR::
SetValue (std::string name, void *src)
{
  // state serialization exceptions
  if (name == "serialization_state") {
    this->load_serialized((char*)src);
    return;
  } else if (name == "serialization_free") {
    this->free_serialized();
    return;
  } else if (name.compare("serialization_create") == 0) {
    this->new_serialized();
    return;
  } else if (name.compare("reset_time") == 0) {
    // time_s and timesteps seems to be used exclusively for reporting current time
    this->state->lgar_bmi_params.time_s = this->GetStartTime();
    this->state->lgar_bmi_params.timesteps = 0;
    return;
  }

  void * dest = NULL;
  dest = this->GetValuePtr(name);

  if (src == NULL) {
    std::stringstream errMsg;
    errMsg << "SetValue: source pointer is null for variable " << name;
    LOG(errMsg.str(), LogLevel::INFO);
    throw std::runtime_error(errMsg.str());
  }

  if (dest == NULL) {
    std::stringstream errMsg;
    errMsg << "SetValue: destination pointer is null for variable " << name;
    LOG(errMsg.str(), LogLevel::INFO);
    throw std::runtime_error(errMsg.str());
  }

  if (name.compare("soil_temperature_profile") == 0) {
    int n = this->state->lgar_bmi_params.num_cells_temp;
    if (n <= 0) {
      std::stringstream errMsg;
      errMsg << "SetValue: invalid num_cells_temp for variable " << name << ": " << n;
      LOG(errMsg.str(), LogLevel::INFO);
      throw std::runtime_error(errMsg.str());
    }

    double *temp = static_cast<double*>(src);
    for (int i = 0; i < n; i++) {
      if (!(temp[i] > 0.0)) {
        std::stringstream errMsg;
        errMsg << "SetValue: soil_temperature_profile[" << i << "] must be > 0.0 K, value=" << temp[i];
        LOG(errMsg.str(), LogLevel::INFO);
        throw std::runtime_error(errMsg.str());
      }
    }
  }

  if (dest) {
    int nbytes = 0;
    nbytes = this->GetVarNbytes(name);
    memcpy(dest, src, nbytes);
  }

}


void BmiLGAR::
SetValueAtIndices (std::string name, int * inds, int len, void *src)
{
  void * dest = NULL;

  dest = this->GetValuePtr(name);

  if (dest) {
    int i;
    int itemsize = 0;
    int offset;
    char *ptr;

    itemsize = this->GetVarItemsize(name);

    for (i=0, ptr=(char *)src; i<len; i++, ptr+=itemsize) {
      offset = inds[i] * itemsize;
      memcpy((char *)dest + offset, ptr, itemsize);
    }
  }
}


std::string BmiLGAR::
GetComponentName()
{
  return "LASAM (Lumped Arid/Semi-arid Model)";
}


int BmiLGAR::
GetInputItemCount()
{
  return this->input_var_name_count;
}


int BmiLGAR::
GetOutputItemCount()
{
  return this->output_var_name_count;
}


std::vector<std::string> BmiLGAR::
GetInputVarNames()
{
  std::vector<std::string> names;

  for (int i=0; i<this->input_var_name_count; i++)
    names.push_back(this->input_var_names[i]);

  return names;
}


std::vector<std::string> BmiLGAR::
GetOutputVarNames()
{
  std::vector<std::string> names;

  for (int i=0; i<this->output_var_name_count; i++)
    names.push_back(this->output_var_names[i]);

  return names;
}


double BmiLGAR::
GetStartTime () {
  return 0.0;
}


double BmiLGAR::
GetEndTime () {
  return this->state->lgar_bmi_params.endtime_s;
}


double BmiLGAR::
GetCurrentTime () {
  return this->state->lgar_bmi_params.time_s;
}


std::string BmiLGAR::
GetTimeUnits() {
  return "s";
}


double BmiLGAR::
GetTimeStep () {
  return this->state->lgar_bmi_params.forcing_resolution_h * 3600.; // convert hours to seconds
}

std::string BmiLGAR::
GetGridType(const int grid)
{
  if (grid == 0)
    return "uniform_rectilinear";
  else
    return "";
}


void BmiLGAR::
GetGridX(const int grid, double *x)
{
  // this is not needed but printing here to avoid compiler warnings
  bmilgar_ss <<"GetGridX: "<<grid<<" "<<x[0]<<"\n";
  LOG(bmilgar_ss.str(), LogLevel::INFO); bmilgar_ss.str("");  
  throw bmi_lgar::NotImplemented();
}


void BmiLGAR::
GetGridY(const int grid, double *y)
{
  // this is not needed but printing here to avoid compiler warnings
  bmilgar_ss <<"GetGridY: "<<grid<<" "<<y[0]<<"\n";
  LOG(bmilgar_ss.str(), LogLevel::INFO); bmilgar_ss.str("");  
  throw bmi_lgar::NotImplemented();
}


void BmiLGAR::
GetGridZ(const int grid, double *z)
{
  // this is not needed but printing here to avoid compiler warnings
  bmilgar_ss <<"GetGridZ: "<<grid<<" "<<z[0]<<"\n";
  LOG(bmilgar_ss.str(), LogLevel::INFO); bmilgar_ss.str("");  
  throw bmi_lgar::NotImplemented();
}


int BmiLGAR::
GetGridNodeCount(const int grid)
{
  // this is not needed but printing here to avoid compiler warnings
  bmilgar_ss <<"GetGridNodeCount: "<<grid<<"\n";
  LOG(bmilgar_ss.str(), LogLevel::INFO); bmilgar_ss.str("");  
  throw bmi_lgar::NotImplemented();
}


int BmiLGAR::
GetGridEdgeCount(const int grid)
{
  // this is not needed but printing here to avoid compiler warnings
  bmilgar_ss <<"GetGridEdgeCount: "<<grid<<"\n";
  LOG(bmilgar_ss.str(), LogLevel::INFO); bmilgar_ss.str("");  
  throw bmi_lgar::NotImplemented();
}


int BmiLGAR::
GetGridFaceCount(const int grid)
{
  // this is not needed but printing here to avoid compiler warnings
  bmilgar_ss <<"GetGridFaceCount: "<<grid<<"\n";
  LOG(bmilgar_ss.str(), LogLevel::INFO); bmilgar_ss.str("");  
  throw bmi_lgar::NotImplemented();
}


void BmiLGAR::
GetGridEdgeNodes(const int grid, int *edge_nodes)
{
  // this is not needed but printing here to avoid compiler warnings
  bmilgar_ss <<"GetGridEdgeNodes: "<<grid<<" "<<edge_nodes[0]<<"\n";
  LOG(bmilgar_ss.str(), LogLevel::INFO); bmilgar_ss.str("");  
  throw bmi_lgar::NotImplemented();
}


void BmiLGAR::
GetGridFaceEdges(const int grid, int *face_edges)
{
  // this is not needed but printing here to avoid compiler warnings
  bmilgar_ss <<"GetGridFaceNodes: "<<grid<<" "<<face_edges[0]<<"\n";
  LOG(bmilgar_ss.str(), LogLevel::INFO); bmilgar_ss.str("");  
  throw bmi_lgar::NotImplemented();
}


void BmiLGAR::
GetGridFaceNodes(const int grid, int *face_nodes)
{
  // this is not needed but printing here to avoid compiler warnings
  bmilgar_ss <<"GetGridFaceNodes: "<<grid<<" "<<face_nodes[0]<<"\n";
  LOG(bmilgar_ss.str(), LogLevel::INFO); bmilgar_ss.str("");  
  throw bmi_lgar::NotImplemented();
}


void BmiLGAR::
GetGridNodesPerFace(const int grid, int *nodes_per_face)
{
  // this is not needed but printing here to avoid compiler warnings
  bmilgar_ss <<"GetGridNodesPerFace: "<<grid<<" "<<nodes_per_face[0]<<"\n";
  LOG(bmilgar_ss.str(), LogLevel::INFO); bmilgar_ss.str("");  
  throw bmi_lgar::NotImplemented();
}


// low-impoact serialization. This will only archive mutable state that is not recalculated during Update or is accessable through BMI getters/setters
template <class Archive>
void BmiLGAR::
serialize(Archive& ar, const unsigned int version) {
  model_state* state = this->state;

  // GetValuePtr properties
  ar & state->lgar_bmi_input_params->precipitation_mm_per_h;
  ar & this->bmi_unit_conv.volprecip_timestep_m;
  ar & state->lgar_bmi_input_params->PET_mm_per_h;
  ar & this->bmi_unit_conv.volPET_timestep_m;
  ar & this->bmi_unit_conv.volAET_timestep_m;
  ar & this->bmi_unit_conv.volrunoff_giuh_timestep_m;
  ar & this->bmi_unit_conv.volend_timestep_m;
  ar & this->bmi_unit_conv.volQ_timestep_m;
  ar & this->bmi_unit_conv.volin_timestep_m;
  ar & this->bmi_unit_conv.volrech_timestep_m;
  ar & this->bmi_unit_conv.volQ_gw_timestep_m;
  ar & this->bmi_unit_conv.mass_balance_m;
  ar & this->bmi_unit_conv.volrunoff_timestep_m;
  ar & this->bmi_unit_conv.volrunoff_giuh_ponded_m;
  ar & state->lgar_bmi_params.num_wetting_fronts;
  ar & state->lgar_calib_params.ponded_depth_max;
  ar & state->lgar_calib_params.field_capacity_psi;

  // end of update
  ar & state->lgar_mass_balance.volprecip_timestep_cm;
  ar & state->lgar_mass_balance.volin_timestep_cm;
  ar & state->lgar_mass_balance.volon_timestep_cm;
  ar & state->lgar_mass_balance.volend_timestep_cm;
  ar & state->lgar_mass_balance.volAET_timestep_cm;
  ar & state->lgar_mass_balance.volrech_timestep_cm;
  ar & state->lgar_mass_balance.volrunoff_timestep_cm;
  ar & state->lgar_mass_balance.volQ_timestep_cm;
  ar & state->lgar_mass_balance.volQ_gw_timestep_cm;
  ar & state->lgar_mass_balance.volPET_timestep_cm;
  ar & state->lgar_mass_balance.volrunoff_giuh_timestep_cm;

  ar & state->lgar_mass_balance.volprecip_cm;
  ar & state->lgar_mass_balance.volin_cm;
  ar & state->lgar_mass_balance.volon_cm;
  ar & state->lgar_mass_balance.volend_cm;
  ar & state->lgar_mass_balance.volAET_cm;
  ar & state->lgar_mass_balance.volrech_cm;
  ar & state->lgar_mass_balance.volrunoff_cm;
  ar & state->lgar_mass_balance.volQ_cm;
  ar & state->lgar_mass_balance.volQ_gw_cm;
  ar & state->lgar_mass_balance.volPET_cm;
  ar & state->lgar_mass_balance.volrunoff_giuh_cm;
  ar & state->lgar_mass_balance.volchange_calib_cm;

  ar & boost::serialization::make_array(state->lgar_bmi_params.cum_layer_thickness_cm, state->lgar_bmi_params.num_layers);

  // giuh state
  ar & boost::serialization::make_array(this->giuh_runoff_queue, state->lgar_bmi_params.num_giuh_ordinates);

  // in frozen_factor_hydraulic_conductivity
  if (state->lgar_bmi_params.sft_coupled) {
    ar & boost::serialization::make_array(state->lgar_bmi_params.frozen_factor, state->lgar_bmi_params.num_layers);
  }

  // update_calibratable_parameters
  ar & state->lgar_bmi_params.field_capacity_psi_cm;
  ar & state->lgar_bmi_params.calib_params_flag;

  // may be set in adapative timesteps
  if (state->lgar_bmi_params.adaptive_timestep){
    ar & state->lgar_bmi_params.timestep_h;
  }

  // how much time has passed since instantiation
  ar & state->lgar_bmi_params.time_s;
  ar & state->lgar_bmi_params.timesteps;

  // serialization of arbitrarily-lengthed linked-lists
  int num_fronts = state->lgar_bmi_params.num_wetting_fronts;
  this->serialize_wetting_front_list(ar, &state->head, state->lgar_bmi_params.num_wetting_fronts);
  int num_previous = 0;
  this->serialize_wetting_front_list(ar, &state->state_previous, num_previous);

  if (Archive::is_loading::value && num_fronts != state->lgar_bmi_params.num_wetting_fronts) {
    // reallocate arrays based on new num_wetting_fronts
    this->realloc_soil();
  }
  ar & boost::serialization::make_array(
    state->lgar_bmi_params.soil_moisture_wetting_fronts, state->lgar_bmi_params.num_wetting_fronts
  );
  ar & boost::serialization::make_array(
    state->lgar_bmi_params.soil_depth_wetting_fronts, state->lgar_bmi_params.num_wetting_fronts
  );

}

template <class Archive>
void BmiLGAR::serialize_wetting_front_list(Archive &ar, wetting_front **head, int &count) {
  wetting_front *current;
  if (Archive::is_saving::value) {
    if (count == 0) // assume recalculation needed if 0
      count = listLength(*head);
    ar & count;
    current = *head;
    while (current != NULL) {
      ar & (*current);
      current = current->next;
    }
  } else { // loading
    ar & count;
    listDelete(*head);
    *head = NULL;
    wetting_front *prior;
    for (int i = 0; i < count; ++i) {
      current = new wetting_front();
      ar & (*current);
      if (i == 0) {
        *head = current;
      } else {
        prior->next = current;
      }
      prior = current;
    }
  }
}

void BmiLGAR::new_serialized() {
  // resize with reserved space for storing size
  this->m_serialized.resize(sizeof(uint64_t));
  boost::archive::binary_oarchive archive(this->m_serialized);
  try {
    archive << (*this);
    this->m_serialized_length = this->m_serialized.size();
    // get serialized size without header and copy size to the beginning of the buffer
    uint64_t serialized_size = this->m_serialized_length - sizeof(uint64_t);
    memcpy(this->m_serialized.data(), &serialized_size, sizeof(uint64_t));
  } catch (const std::exception &e) {
    LOG(LogLevel::INFO, "Serializing LASAM encountered an error: %s", e.what());
    this->free_serialized();
    throw;
  }
}

void BmiLGAR::load_serialized(char* data) {
  // copy size from the start of data
  uint64_t size;
  memcpy(&size, data, sizeof(uint64_t));
  // create stream from everything past the size header
  membuf stream(data + sizeof(uint64_t), size);
  boost::archive::binary_iarchive archive(stream);
  try {
    archive >> (*this);
  } catch (const std::exception &e) {
    LOG(LogLevel::INFO, "Deserializing LASAM encountered an error: %s", e.what());
    throw;
  }
  this->free_serialized();
}

void BmiLGAR::free_serialized() {
  this->m_serialized.clear();
  this->m_serialized.shrink_to_fit();
  this->m_serialized_length = 0;
}


#endif
