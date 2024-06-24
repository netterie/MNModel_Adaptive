# ****************************************************************************
# Code to understand the specification of parameters
# ****************************************************************************

# ----------------------------------------------------------------------------
# Setup
# ----------------------------------------------------------------------------

# Set directories based on root set by "COVID19_Trigger.RProj"
dir_root  <- getwd()
dir_proj  <- file.path(dir_root, "COVID19Proj")
dir_data  <- file.path(dir_proj, "data")
dir_setup <- file.path(dir_proj, "setup")

# Empty workspace and load libraries 
library(COVID19pack) # make sure you have installed the COVID19pack package
library(ggplot2)
library(reshape2)

# ----------------------------------------------------------------------------
# gen_abc_param23 (function from COVID19pack/R/spec_abc.R)
# ----------------------------------------------------------------------------

# In Scenarios.R, the "parms" object is initialized by calling "gen_abc_param23".
# Zongbo explained that there were many iterations of calibration and 
# this was the set of calibration parameters that they finally settled on


# INITIALIZE -----------------------------------------------------------------

# From Scenarios.R: 
# ---Source helper functions in the setup folder
setup_files = as.character(file.path(dir_setup, list.files(dir_setup)))
sapply(setup_files, source)
# ---Load the calibrated samples and use j=1 (first set) to explore
x<-readRDS(file.path(dir_data, "abc_posterior_samples.rds"))
j<-1
# ---Load ..defaults? into environment
tmp_ls <- set_states_eo(as.Date("2021-03-22"))
list2env(tmp_ls, envir = .GlobalEnv)

# Initialize inputs to "gen_abc_param23" according to the following:
# parms <- COVID19pack:::gen_abc_param23(x[j,1:17], 
#                                        epi_states_index = epi_states_index, 
#                                        hash_table = hash_table)
# gen_abc_param23 <- function(x, epi_states_index, hash_table,
#                             mask_cr_mat = NULL, # estimated contact reduction matrix including masking
#                             start_time_mask_policy = Inf,
#                             beta_after = 0.02, 
#                             beta_change = FALSE,
#                             hosp_change = FALSE,
#                             hosp_mul = 1) {

x <- x[j,1:17]
#epi_states_index = epi_states_index
#hash_table = hash_table
mask_cr_mat = NULL
start_time_mask_policy = Inf
beta_after = 0.02
beta_change = FALSE
hosp_change = FALSE
hosp_mul = 1

# RUN FUNCTION CODE -------------------------------------------------------------

  sip_cr_mat <- matrix(c(x[[9]],  x[[10]], x[[10]], x[[11]],
                         x[[10]], x[[12]], x[[12]], x[[11]],
                         x[[10]], x[[12]], x[[12]], x[[11]],
                         x[[11]], x[[11]], x[[11]], x[[13]]),
                       nrow = 4, byrow = T)
  
  beh_cr_mat <- matrix(c(x[[14]]*x[[9]], x[[15]],         x[[15]],         x[[17]],
                         x[[15]],        x[[16]]*x[[12]], x[[16]]*x[[12]], x[[17]],
                         x[[15]],        x[[16]]*x[[12]], x[[16]]*x[[12]], x[[17]],
                         x[[17]],        x[[17]],         x[[17]],         x[[13]]),
                       nrow = 4, byrow = T)

# ----------------------------------------------------------------------------
# parameters.R (function from COVID19pack/R/spec_abc.R)
# ----------------------------------------------------------------------------

# INITIALIZE -----------------------------------------------------------------
# (this is done within gen_abc_param23)

  # Initialize the function inputs NOT specified in gen_abc_param23 that get 
  # set to their defaults. I am listing all inputs and just commenting out
  # those that are defined above via gen_abc_param23 inputs
  
  #beta = 0.019,       # from calibration
  #beta_after = 0.019, # from calibration
  #beta_change = FALSE, 
  #hosp_mul = 1, 
  #hosp_change = FALSE, 
  #weight_60p_init = 0.08 # from calibration
  n_days_incubation = 5.2
  n_days_infectious = 7.8
  exposed_transition_rate = NA
  infected_transition_rate = NA
  n_days_rec_hosp = c(4.41, 4.41, 4.41, 5.80, 5.80, 
                      6.03, 7.72, 7.60, 6.45) # days by age group
  n_days_rec_ICU = c(17.50, 17.50, 17.50, 17.50, 17.50, 
                     19.65, 21.51, 15.38, 15.38) # days by age group                       
  #init_cases_detected =  0.019, # from calibration
  n_icu_beds = 2200
  relative_risk_mort_co = 1 # Same risk with comorbidity (e.g. comorbidities not modeled)
  str_peak_type = c("hospitalizations") # or "deaths" or "infections" or "new_cases", "new_sym_cases"
  v_strat_status = c("sd" = 1, "sip" = 1, "sc" = 1, "sd60p" = 1, "bec" = 1)
  prop_asymptomatic = NULL # if NULL, we assume there is an age-difference for asymptomatic infection
  #prop_asymp_19 = 0.408, # from calibration, proportion 0-19 year-old asymptomatic
  #p_h_50_59 = c(0.062, 0.062), # from calibration
  #p_h_60p = c(0.245, 0.245), # from calibration
  p_dying_home_70 = 0.088 # from calibration
  #p_dying_home_80 = 0.088, # from calibration
  ##### Parameters related to social distancing ####
  start_time_social_distancing = Inf # set when running the model
  start_time_sip = Inf
  start_time_behavior_change = Inf
  start_time_60plus_distancing = Inf # NOT USED
  start_time_school_closure = Inf    # NOT USED
  end_time_social_distancing = -Inf # set when running the model
  end_time_sip = -Inf
  end_time_behavior_change = -Inf
  end_time_60plus_distancing = -Inf # NOT USED
  end_time_school_closure = -Inf    # NOT USED
  #social_distancing_contact_reduction = 0.246,   # from calibration, contract reduction in period 1
  #sip_cr_mat = matrix(c(0.65, 0.65, 0.65, 0.148, # from calibration, contact reduction in period 2
  #                      0.65, 0.65, 0.65, 0.148, 
  #                      0.65, 0.65, 0.65, 0.148, 
  #                      0.148, 0.148, 0.148, 0.062), 
  #                    nrow = 4, byrow = T),
  #beh_cr_mat = matrix(c(0.454, 0.931, 0.454, 0.735, 
  #                      0.931, 0.291, 0.931, 0.931, 
  #                      0.454, 0.931, 0.454, 0.735, 
  #                      0.735, 0.931, 0.735, 0.062), 
  #                    nrow = 4, byrow = T), # from calibration, contact reduction in period 3
  school_closure_contact_reduction = 0.4 # not used
  sip_days_past_peak = -Inf
  social_distancing_days_past_peak = -Inf
  sixty_plus_days_past_peak = -Inf
  school_closure_days_past_peak = -Inf
  behavior_change_days_past_peak = -Inf
  #### Parameters related to testing (NOT USED) ####
  p_sens_test = 0.8
  p_spec_test = 0.99
  n_tests_per_day = 0
  quarantine_contact_reduction = 0.8
  n_days_quarantine = 14
  p_I_seek_test = 0.05
  p_nonI_seek_test = 0.0025
  n_tests_per_day_hosp = 2/11
  n_tests_per_day_icu = 2/8
  start_time_testing = 56
  ##### Cycling, or trigger, strategies ####
  trigger_strategy = NULL 
  ##### Parameters related to hospitalizations post May 15, 2020 ("Medical Advances") ####
  p_decrease_hosp_mort = 0 # no longer multiplier; estimated from data
  p_decrease_icu_mort = 0  # no longer multiplier; estimated from data
  n_days_rec_hosp_ma = c(3.85, 3.85, 3.85, 5.02, 5.02, 
                         6.45, 7.96, 6.64, 6.55) # by age group
  n_days_rec_ICU_ma = c(15.33, 15.33, 15.33, 15.33, 15.33, 
                        16.86, 19.10, 15.47, 15.47) # by age group
  start_time_ma_mort = Inf # Mortality held constant over simulation period (no impact of "Medical Advances")
  start_time_ma_dur = as.numeric(as.Date("2020-05-15") - set_start_date()) # Change in length of stay post May 15, 2020 ("Medical Advances")
  ##### Masking (NOT USED) ####
  mask_cr_mat = NULL # estimated contact reduction matrix including masking
  start_time_mask_policy = Inf # the starting time of the masking policy 
  ##### Timestep ####
  timestep = 0.05 # note, time step must be 1 or less
  
# Initialize the function inputs specified within gen_abc_param23
#  parameters(epi_states_index = epi_states_index,
#             hash_table = hash_table,
             init_cases_detected = x[[1]]
             p_h_50_59 = c(x[[2]], x[[2]])
             p_h_60p = c(x[[3]], x[[3]])
             p_dying_home_80 = x[[4]]
             prop_asymp_19 = x[[5]]
             beta = x[[6]]
             weight_60p_init = x[[7]]
             social_distancing_contact_reduction = x[[8]]
             sip_cr_mat = sip_cr_mat
             beh_cr_mat = beh_cr_mat
             mask_cr_mat = mask_cr_mat #NULL
             start_time_mask_policy = start_time_mask_policy #Inf
             beta_after = beta_after #0.02
             beta_change = beta_change #FALSE
             hosp_mul = hosp_mul #1
             hosp_change = hosp_change #FALSE
             

             
# RUN FUNCTION CODE -------------------------------------------------------------
          
   
  # Changing the "for (x...)" to "for (i...)" so as to not get mixed up with 
  # the x that is the calibrated data
  # epi_states_index is a list of model structure parameters and also the master
  # index vector "ls_ix". This code loads those list elements into the environment
  for (i in c(1:length(epi_states_index))) {
    assign(names(epi_states_index)[i], epi_states_index[[i]])
  }
  
  # This is TRUE and this code does run
  # It updates the prop_asymptomatic MATRIX based on the calibrated value of prop_asymp_19
  # The matrix has two rows for each comorbidity group and 9 columns for each age group
  # The code sets both rows to the same values of the vector "tmp"
  if (is.null(prop_asymptomatic)) {
    prop_asymptomatic <- matrix(0, nrow = ncg, ncol = n_age_groups)
    tmp <- c(1, 1.113, 1.028, 0.944, 0.845, 0.718, 0.521, 0.437, 0.437) * prop_asymp_19  # Davies et al.
    prop_asymptomatic[1, ] <- prop_asymptomatic[2, ] <- tmp
  } else {
    prop_asymptomatic <- matrix(prop_asymptomatic, nrow = 2, ncol = 9)
  }

  ## proportion of people hospitalized by age and time period
  v_p_h1 <- c(rep(p_h_50_59[1], 6), rep(p_h_60p[1], 3))
  v_p_h2 <- c(rep(p_h_50_59[2], 6), rep(p_h_60p[2], 3))
  prop_hosp_ls <- list(rep(v_p_h1 * c(0.16, 0.16, 0.27, 0.40, 0.55, 1, 1, 1, 1), 
                           each = ncg), # 3/22/2020 - 5/31/2020
                       rep(v_p_h2 * c(0.26, 0.26, 0.16, 0.32, 0.63, 1, 1, 1, 1), 
                           each = ncg)) # 6/1/2020 - present
  
  ## proportion of hospitalized people requiring mechanical ventilation by age and time period
  prop_ICU_ls <- list(rep(c(0.0219, 0.0219, 0.0219, 0.2069, 0.1864,
                            0.3053, 0.3333, 0.2198, 0.1429), each = ncg), # 3/22/2020-4/14/2020
                      rep(c(0.0219, 0.0219, 0.0219, 0.0897, 0.1529, 
                            0.1672, 0.1344, 0.1399, 0.0656), each = ncg), # 4/15/2020-5/14/2020
                      rep(c(0.0260, 0.0260, 0.0260, 0.0387, 0.0808, 
                            0.0981, 0.1166, 0.1005, 0.0425), each = ncg), # 5/15/2020-6/14/2020
                      rep(c(0.0146, 0.0146, 0.0146, 0.0314, 0.0352, 
                            0.0584, 0.0644, 0.0417, 0.0425), each = ncg)) # 6/15/2020-present
  
  # Proportion with comorbidities by age (effectively NOT USED)
  comorbidity_prop_by_age <- c(1.2, 2.3, 3.8, 6.9, 12.5, 21.4, 30.1, 43.6, 55.9) / 100
  
  ## Exposed and Infected compartment dynamics
  # exposed compartments
  # This is TRUE and the code runs
  if (is.na(exposed_transition_rate)) {
    exposed_transition_rate <- n_exposed_states / n_days_incubation 
  }
  # infected compartments
  if (is.na(infected_transition_rate)) {
    infected_transition_rate <- n_infected_states / n_days_infectious
  }
  
  # Expand parameter matrices which is age or comorbidity dependent (no difference by comorbidity status)
  
  ## n_days_rec_hosp is the number of days a person will spend in the hospital if admitted
  n_days_rec_hosp <- matrix(rep(n_days_rec_hosp, ncg), byrow = TRUE, nrow = ncg)
  n_days_rec_hosp_ma <- matrix(rep(n_days_rec_hosp_ma, ncg), byrow = TRUE, nrow = ncg)
  
  #n_days_rec_ICU is the number of days someone will spend in the ICU if they have a bed
  n_days_rec_ICU <- matrix(rep(n_days_rec_ICU, ncg), byrow = TRUE, nrow = ncg)
  n_days_rec_ICU_ma <- matrix(rep(n_days_rec_ICU_ma, ncg), byrow = TRUE, nrow = ncg)
  
  #icu_death_rate_1 is the rate at which people in the ICU die if they have a bed
  icu_death_rate_1 <- matrix(0, nrow = ncg, ncol = n_age_groups)
  
  ## Probability of dying if requiring mechanical ventilation, by age
  prob_icu_death <- c(0,
                      0,
                      0,
                      0.16,
                      0.1739,
                      0.3053,
                      0.4894,
                      0.6029,
                      0.7632)
  
  icu_death_rate_1[1,] <- -log(1-prob_icu_death)
  
  icu_death_rate_1[2,] <- icu_death_rate_1[1,] * relative_risk_mort_co
  prob_icu_death <- 1 - exp(-1 * icu_death_rate_1)
  prob_icu_death <- prob_icu_death
  
  #assume that people who need an icu bed but don't get one, will always die
  prob_icu_death_no_bed <- matrix(1,nrow = 2, ncol = 9 )
  
  #prop_hosp_die is the proportion of hospitalized people who are hospitalized and die
  prob_hosp_death <- matrix(0, nrow = 2, ncol = 9)
  h_death_rate <- matrix(0, nrow = 2, ncol = 9)
  
  ## Probability of dying if hospitalized (without mechanical ventilation) by age
  prob_hosp_death[1, ] <- c(0,
                            0,
                            0,
                            0.0077,
                            0.0097,
                            0.0526,
                            0.1091,
                            0.2154,
                            0.3974)
  
  prob_hosp_death[2,] <- prob_hosp_death[1,]
  h_death_rate <- -log(1 - prob_hosp_death)
  h_death_rate[2,] <- h_death_rate[1,] * relative_risk_mort_co
  prob_hosp_death <- 1 - exp(-1 * h_death_rate)
  
  # Probability of symptomatic infectious dying out of hospital (only oldest age groups)
  prop_inf_die <- matrix(0, nrow = 2, ncol = 9)
  prop_inf_die[, 8] <- p_dying_home_80
  prop_inf_die[, 9] <- p_dying_home_80
  prop_inf_die[2,8:9] <- -log(1 - prop_inf_die[2,8:9])*relative_risk_mort_co
  prop_inf_die[2,8:9] <- 1 - exp(-1 * prop_inf_die[2,8:9])
  
  # Proportion of population in each age group
  prop_inf_by_age <- c(0.08930854, 0.26353432, 0.15794502, 0.15717893, 0.14297132, 0.11012012, 0.04463342, 0.02264691, 0.01166143)
  # Readjust infections for proportion infections in 60+ ages
  if (!is.null(weight_60p_init)) {
    prop_inf_by_age <- c((1 - weight_60p_init) * prop_inf_by_age[1:6] / sum(prop_inf_by_age[1:6]), 
                         weight_60p_init * prop_inf_by_age[7:9] / sum(prop_inf_by_age[7:9]))
  }
  
  #### Contact reduction matrix under SAH: 
  #### Turn a 4x4 indicating the values for 4 age blocks into 
  # Turn it into a vector
  sip_flat <- c(t(sip_cr_mat)) 
  # The "m" matrices are 0/1 age x age indicator matrices that turn
  # on contact reductions for different sets of ages. So slip_flat[x] 
  # is a scalar, and hash_table$m[[x]] is a matrix. 
  # sip_contact_reduction ends up being a list of matrices that have
  # zeros for most of the entries except for where the sip_flat scalar
  # was multiplied by 
  sip_contact_reduction <- lapply(c(1:length(sip_flat)), function(x) {
    sip_flat[x] * hash_table$m[[x]]
  })
  # End result is an age x age matrix that is the sum of all the
  # matrices created in the previous line.
  # JKB Note: I think it's better just to specify these directly. This
  # is a lot of difficult-to-follow code just to accomplish something
  # that could be quickly hand-coded in Excel
  sip_contact_reduction <- Reduce("+", sip_contact_reduction)
  
  #### Contact reduction post SAH
  beh_cr_flat <- c(t(beh_cr_mat))
  behavior_change_contact_reduction <- lapply(c(1:length(sip_flat)), function(x) {
    beh_cr_flat[x] * hash_table$m[[x]]
  })
  behavior_change_contact_reduction_o <- Reduce("+", behavior_change_contact_reduction)
  
  behavior_change_contact_reduction <- behavior_change_contact_reduction_o
  behavior_change_contact_reduction[7:9, ] <- 0
  behavior_change_contact_reduction[, 7:9] <- 0
  
  #### Adjust contact rate of each age group based on what social distancing practices are in place
  school_closure_contact_reduction <- school_closure_contact_reduction * hash_table$m_school
  sixty_plus_contact_reduction <- behavior_change_contact_reduction_o
  sixty_plus_contact_reduction[1:6, 1:6] <- 0
  
  #### Don't adjust Beta (if applicable)
  if (beta_change == FALSE) {
    beta_after <- beta
  }
  #### Don't adjust hospitalization risk (if applicable)
  if (hosp_change == FALSE) {
    hosp_mul <- 1
  }
  
  
  parms <- list("beta" = beta, 
                "beta_after" = beta_after, 
                "hosp_mul" = hosp_mul,
                "n_exposed_states" = n_exposed_states,
                "n_infected_states" = n_infected_states,
                "init_cases_detected" = init_cases_detected,
                "N" = N,
                "N_by_age" = N_by_age,
                "n_icu_beds" = n_icu_beds,
                "n_age_groups" = n_age_groups,
                "n_epi_groups" = n_epi_groups,
                "n_co_groups" = ncg,
                "start_time_social_distancing" = start_time_social_distancing,
                "start_time_school_closure" = start_time_school_closure,
                "start_time_sip" = start_time_sip,
                "start_time_60plus_distancing"=start_time_60plus_distancing,
                "start_time_behavior_change" = start_time_behavior_change,
                "end_time_60plus_distancing"=end_time_60plus_distancing,
                "end_time_social_distancing" = end_time_social_distancing,
                "end_time_school_closure" = end_time_school_closure,
                "end_time_sip" = end_time_sip,
                "end_time_behavior_change" = end_time_behavior_change,
                "sixty_plus_contact_reduction" = sixty_plus_contact_reduction,
                "sip_days_past_peak" = sip_days_past_peak,
                "social_distancing_days_past_peak" = social_distancing_days_past_peak,
                "sixty_plus_days_past_peak" = sixty_plus_days_past_peak,
                "school_closure_days_past_peak" = school_closure_days_past_peak,
                "behavior_change_days_past_peak" = behavior_change_days_past_peak,
                "str_peak_type" = str_peak_type,
                "age_groups" = age_groups,
                "timestep" = timestep,
                "sip_contact_reduction" = sip_contact_reduction,
                "behavior_change_contact_reduction_o" = behavior_change_contact_reduction_o, 
                "behavior_change_contact_reduction" = behavior_change_contact_reduction, 
                "social_distancing_contact_reduction" = social_distancing_contact_reduction,
                "school_closure_contact_reduction" = school_closure_contact_reduction,
                "age_prop" = age_prop,
                "comorbidity_prop_by_age" = comorbidity_prop_by_age,
                "v_strat_status" = v_strat_status,
                "ls_ix" = ls_ix, # index list
                "epi_groups_ls" = epi_groups_ls, # replace the separate vector states (v_exp_str, v_qe_str, v_inf_str, v_asym_inf_str, v_qi_str, and v_qai_str)
                "prop_asymptomatic" = prop_asymptomatic, # flatten the matrix: the matrix is flattened by columns
                "prop_inf_die" = prop_inf_die,
                "prop_inf_by_age" = prop_inf_by_age,
                "prob_icu_death_no_bed" = prob_icu_death_no_bed,
                "p_sens_test" = p_sens_test,
                "p_spec_test" = p_spec_test,
                "hash_table" = hash_table,
                "quarantine_contact_reduction" = quarantine_contact_reduction,
                "start_time_testing" = start_time_testing,
                "p_decrease_hosp_mort" = p_decrease_hosp_mort, 
                "p_decrease_icu_mort" = p_decrease_icu_mort,
                "start_time_ma_mort" = start_time_ma_mort, 
                "start_time_ma_dur" = start_time_ma_dur, 
                #####cycling and trigger strategies####
                "trigger_strategy" = trigger_strategy, 
                ##### Masking ########
                "mask_cr_mat" = mask_cr_mat, 
                "start_time_mask_policy" = start_time_mask_policy, 
                #####whether reduction parameter is modified####
                "reduction_modified" = 0)
  
  init_vec <- get_initial_conditions(parms, m_init_cases)
  names(init_vec) <- nam
  parms$init_vec <- init_vec
  
  ############################################################
  #### Flatten matrices for difference vector calculation #### 
  ############################################################
  
  # JKB: This part is actually just initializing VECTOR parameters that
  # change over time to the iniial value, and setting the 
  # date that they change. 
  parms$prop_hosp_ls <- prop_hosp_ls
  parms$prop_hosp <- prop_hosp_ls[[1]]
  parms$prop_hosp_change_day <- as.numeric(c(set_start_date(), 
                                            as.Date(c("2020-06-01"))) - 
                                            set_start_date())
  
  parms$prop_ICU_ls <- prop_ICU_ls
  parms$prop_ICU <- prop_ICU_ls[[1]] 
  parms$prop_ICU_change_day <- as.numeric(c(set_start_date(), 
                                            as.Date(c("2020-04-15", "2020-05-15", "2020-06-15"))) - 
                                            set_start_date())
  
  # This part is flatting matrices that had a row for each comorbidity. 
  # Now they are vectors where age varies fastest and comorbidity varies slower
  parms$prob_icu_death_no_bed <- c(prob_icu_death_no_bed)
  parms$prop_inf_die <- c(prop_inf_die)
  parms$prop_asymptomatic <- c(prop_asymptomatic)
  
  ###########################################
  #### Converting and rescale parameters ####
  ###########################################
  ## This part includes code in convert_parms.R and 
  ## parameter calculation in covid_19_model_function_v3.R
  
  timestep <- parms$timestep
  
  # JKB: Eva suggested that specifying the parameters
  # correctly per time step should be a requirement of
  # the input data, rather than being integrated into the code
  n_days_rec_ICU <- n_days_rec_ICU / timestep
  n_days_rec_ICU_ma <- n_days_rec_ICU_ma / timestep
  n_days_rec_hosp <- n_days_rec_hosp / timestep
  n_days_rec_hosp_ma <- n_days_rec_hosp_ma / timestep
  n_days_quarantine <- n_days_quarantine / timestep
  
  parms$p_icu_exit <- list("null" = c(1 - exp(-1 / n_days_rec_ICU)), 
                           "MedAdv" = c(1 - exp(-1 / (n_days_rec_ICU_ma))))
  parms$p_h_exit <- list("null" = c(1 - exp(-1 / n_days_rec_hosp)), 
                         "MedAdv" = c(1 - exp(-1 / (n_days_rec_hosp_ma))))
  parms$p_quar_exit <- 1 - exp(-1 / n_days_quarantine)
  parms$p_icu_exit_nobed <- 1 - exp( - timestep )
  
  parms$prob_hosp_death <- list("null" = c(prob_hosp_death), 
                                "MedAdv" = c(prob_hosp_death * (1 - p_decrease_hosp_mort)))
  parms$prob_icu_death <- list("null" = c(prob_icu_death), 
                               "MedAdv" = c(prob_icu_death * (1 - p_decrease_icu_mort)))
  
  exposed_transition_rate <- exposed_transition_rate * timestep
  infected_transition_rate <- infected_transition_rate * timestep
  parms$p_trans_exp <- 1 - exp(-1 * exposed_transition_rate)
  parms$p_trans_inf <- 1 - exp(-1 * infected_transition_rate)
  parms$exposed_transition_rate <- exposed_transition_rate
  parms$infected_transition_rate <- infected_transition_rate
  
  # mixing_matrix_o is the original mixing matrix (unmodified)
  # mixing_matrix is the matrix that would be modified if a peak or a date is detected
  parms$mixing_matrix_o <- parms$mixing_matrix <- mixing_matrix * timestep
  
  parms$n_tests_per_day <- n_tests_per_day * timestep
  parms$n_tests_per_day_hosp <- n_tests_per_day_hosp * timestep
  parms$n_tests_per_day_icu <- n_tests_per_day_icu * timestep
  
  parms$p_I_seek_test <- 1 - exp(log(1 - p_I_seek_test) * timestep)
  parms$p_nonI_seek_test <- 1 - exp(log(1 - p_nonI_seek_test) * timestep)
