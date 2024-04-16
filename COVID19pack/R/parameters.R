#' @title Model parameter generation
#'
#' Parameter Generation for MN Covid-19 model
#'
#' input values of parameters or use defaults to generate list of parameters
#' that is used by main function.
#' 
#' @param epi_states_index a list with sub-epi states and indices. This list is generated from `gen_states_index()`. 
#' @param hash_table a hash table for indexing different epistates, age group and comorbidity group throught a simulation.
#' @param p_sens_test sensitivity for test-probability of people who test positive among people infected.
#' @param p_spec_test specificity for test-probability of people who test negative among people not infected.
#' @param n_tests_per_day number of tests per day.
#' @param quarantine_contact_reduction the percentage reduction in contact rates when quarantine in effect.
#' @param n_days_quarantine number of days that quarantine in effect.
#' @param p_I_seek_test probability of infected individuals seeking test.
#' @param p_nonI_seek_test  probability of non-infected individuals seeking test. 
#' @param n_tests_per_day_hosp the number of hospitalized individuals being tested per day.
#' @param n_tests_per_day_icu the number of individuals in ICU being tested per day.
#' @param start_time_testing day of covid-19 testing start. 
#' @param trigger_strategy specify the trigger_strategy. The default is `NULL`
#' @param p_decrease_hosp_mort percentage of decrease in probability of dying at hospital.
#' @param p_decrease_icu_mort  percentage of decrease in probability of dying at ICU.
#' @param p_decrease_hosp_dur duration of decrease in probability of dying at hospital.
#' @param p_decrease_icu_dur duration of decrease in probability of dying at ICU.
#' @param start_time_ma_mort day of medical advances start to reduce mortality.
#' @param start_time_ma_dur day of medical advances start to reduce duration in hospital and ICU bed.
#' @param age_groups vector which lists the starting age of each age group.
#' @param beta probability that contact with an infectious person will cause an infection.
#' @param n_days_incubation average number of days spent in the exposed states.
#' @param n_days_infectious average number of days spent in the infectious states.
#' @param exposed_transition_rate rate at which population moves through the exposed tunnel states. This value is NA by default,
#'  and if it is not supplied directly, the value is calculated internally in accordance with timestep, n_exposed_states, and n_days_incubation.
#' @param infected_transition_rate rate at which population moves through the exposed tunnel states. This value is NA by default,
#'  and if it is not supplied directly, the value is calculated internally in accordance with timestep, n_infected_states, and n_days_infectious.
#' @param n_days_rec_hosp average number of days it takes to recover prior to 5/15 if in the hospital.
#' @param n_days_rec_hosp_ma average number of days it takes to recover post 5/15 if in the hospital.
#' @param n_days_rec_ICU average number of days it takes to recover if in the ICU (w/ bed) prior to 5/15.
#' @param n_days_rec_ICU_ma average number of days it takes to recover if in the ICU (w/ bed) post 5/15 if in the ICU.
#' @param init_cases_detected percentage of cases that are assumed to have been detected, as of day 0 in the model. 
#' Used in inital condition calculation.
#' @param n_icu_beds number of ICU beds available in Minnesota.
#' @param relative_risk_mort_co relative risk of mortality if a person has 1 or more comorbidities.
#' @param start_time_social_distancing day of simulation at which social distancing measures start.
#' @param start_time_school_closure day of simulation at which school closures start.
#' @param start_time_sip day of simulation at which shelter in place starts.
#' @param start_time_60plus_distancing day of simulation at which targeted social distancing measures start for individuals age 60+.
#' @param start_time_behavior_change day of simulation at which population changes behavior independent of enacted policies to start socially distancing.
#' @param end_time_60plus_distancing day of simulation at which targeted social distancing measures end for individuals age 60+.
#' @param end_time_social_distancing day of simulation at which social distancing measures end.
#' @param end_time_school_closure day of simulation at which school closures end.
#' @param end_time_sip day of simulation at which shelter in place ends.
#' @param end_time_behavior_change day of simulation at which population ceases policy-independent social distancing.
#' @param sip_cr_mat as a 4x4 matrix with 4 age groups: 0-19, 20-39, 40-59, 60+. This matrix represents the contact reduction parameter between two age groups or within an age group in SAH. 
#' @param beh_cr_mat as a 4x4 matrix with 4 age groups: 0-19, 20-39, 40-59, 60+. This matrix represents the contact reduction parameter between two age groups or within an age group post SAH. 
#' @param social_distancing_contact_reduction as a decimal, the percentage reduction in contact rates when social distancing is in effect
#' @param school_closure_contact_reduction as a decimal, the percentage reduction in contact rates when school closures are in effect. 
#' This reduction only applies to relevant entries of the mixing matrix.
#' @param sip_days_past_peak number of days shelter in place will be in effect past the chosen peak (infections, hospitalizations, deaths)
#' @param social_distancing_days_past_peak number of days social distancing will be in effect past the chosen peak (infections, hospitalizations, deaths)
#' @param sixty_plus_days_past_peak number of days 60+ social distancing will be in effect past the chosen peak (infections, hospitalizations, deaths)
#' @param school_closure_days_past_peak number of days schools remain closed past the chosen peak (infections, hospitalizations, deaths)
#' @param behavior_change_days_past_peak number of days a long-lasting behavior change will be in effect past the chosen peak (infections, hospitalizations, deaths)
#' @param str_peak_type string indicating which peak type will be used in calculations. The string could be "deaths", "hospitalizations", "infections", "new_cases", or "new_sym_cases". 
#' @param v_strat_status a vector that will be used in the solve_model() function which indicates which social distancing strategies are being used
#' @param prop_asymptomatic vector that contains the proportion of each age groups which is asymptomatic
#' @param p_h_50_59 a vector of two elements, representing the probability of hospitalization among age 50-59. The first element is the probability in April/May and the second element is the probability in June-present. 
#' @param p_h_60p a vector of two elements, representing the probability of hospitalization among age 60+. The first element is the probability in April/May and the second element is the probability in June-present. 
#' @param p_dying_home_70 probability of dying at home for 70-79 year old people who have a symptomatic infection
#' @param p_dying_home_80 probability of dying at home for 80+ year old people who have a symptomatic infection
#' @param mask_cr_mat a matrix of contact reduction considering/including masking effect. The default is `NULL`. 
#' @param start_time_mask_policy The start time of masking policy. The default is `Inf`. 
#' @param timestep length of timesteps in days used by the solver function, should be 1 or less
#'
#' @return
#' list of parameters
#'
#' @export
parameters <- function(epi_states_index, 
                       hash_table, 
                       beta = 0.019,       # from calibration
                       beta_after = 0.019, # from calibration
                       beta_change = FALSE, 
                       hosp_mul = 1, 
                       hosp_change = FALSE, 
                       weight_60p_init = 0.08, # from calibration
                       n_days_incubation = 5.2,
                       n_days_infectious = 7.8,
                       exposed_transition_rate = NA,
                       infected_transition_rate = NA,
                       n_days_rec_hosp = c(4.41, 4.41, 4.41, 5.80, 5.80, 
                                           6.03, 7.72, 7.60, 6.45), # days by age group
                       n_days_rec_ICU = c(17.50, 17.50, 17.50, 17.50, 17.50, 
                                          19.65, 21.51, 15.38, 15.38), # days by age group                       
                       init_cases_detected =  0.019, # from calibration
                       n_icu_beds = 2200,
                       relative_risk_mort_co = 1, # Same risk with comorbidity (e.g. comorbidities not modeled)
                       str_peak_type = c("hospitalizations"), # or "deaths" or "infections" or "new_cases", "new_sym_cases"
                       v_strat_status = c("sd" = 1, "sip" = 1, "sc" = 1, "sd60p" = 1, "bec" = 1),
                       prop_asymptomatic = NULL, # if NULL, we assume there is an age-difference for asymptomatic infection
                       prop_asymp_19 = 0.408, # from calibration, proportion 0-19 year-old asymptomatic
                       p_h_50_59 = c(0.062, 0.062), # from calibration
                       p_h_60p = c(0.245, 0.245), # from calibration
                       p_dying_home_70 = 0.088, # from calibration
                       p_dying_home_80 = 0.088, # from calibration
                       ##### Parameters related to social distancing ####
                       start_time_social_distancing = Inf, # set when running the model
                       start_time_sip = Inf,
                       start_time_behavior_change = Inf,
                       start_time_60plus_distancing = Inf, # NOT USED
                       start_time_school_closure = Inf,    # NOT USED
                       end_time_social_distancing = -Inf, # set when running the model
                       end_time_sip = -Inf,
                       end_time_behavior_change = -Inf,
                       end_time_60plus_distancing = -Inf, # NOT USED
                       end_time_school_closure = -Inf,    # NOT USED
                       social_distancing_contact_reduction = 0.246,   # from calibration, contract reduction in period 1
                       sip_cr_mat = matrix(c(0.65, 0.65, 0.65, 0.148, # from calibration, contact reduction in period 2
                                             0.65, 0.65, 0.65, 0.148, 
                                             0.65, 0.65, 0.65, 0.148, 
                                             0.148, 0.148, 0.148, 0.062), 
                                           nrow = 4, byrow = T),
                       beh_cr_mat = matrix(c(0.454, 0.931, 0.454, 0.735, 
                                             0.931, 0.291, 0.931, 0.931, 
                                             0.454, 0.931, 0.454, 0.735, 
                                             0.735, 0.931, 0.735, 0.062), 
                                           nrow = 4, byrow = T), # from calibration, contact reduction in period 3
                       school_closure_contact_reduction = 0.4, # not used
                       sip_days_past_peak = -Inf,
                       social_distancing_days_past_peak = -Inf,
                       sixty_plus_days_past_peak = -Inf,
                       school_closure_days_past_peak = -Inf,
                       behavior_change_days_past_peak = -Inf,
                       #### Parameters related to testing (NOT USED) ####
                       p_sens_test = 0.8,
                       p_spec_test = 0.99,
                       n_tests_per_day = 0,
                       quarantine_contact_reduction = 0.8,
                       n_days_quarantine = 14,
                       p_I_seek_test = 0.05,
                       p_nonI_seek_test = 0.0025,
                       n_tests_per_day_hosp = 2/11,
                       n_tests_per_day_icu = 2/8,
                       start_time_testing = 56,
                       ##### Cycling, or trigger, strategies ####
                       trigger_strategy = NULL, 
                       ##### Parameters related to hospitalizations post May 15, 2020 ("Medical Advances") ####
                       p_decrease_hosp_mort = 0, # no longer multiplier; estimated from data
                       p_decrease_icu_mort = 0,  # no longer multiplier; estimated from data
                       n_days_rec_hosp_ma = c(3.85, 3.85, 3.85, 5.02, 5.02, 
                                              6.45, 7.96, 6.64, 6.55), # by age group
                       n_days_rec_ICU_ma = c(15.33, 15.33, 15.33, 15.33, 15.33, 
                                             16.86, 19.10, 15.47, 15.47), # by age group
                       start_time_ma_mort = Inf, # Mortality held constant over simulation period (no impact of "Medical Advances")
                       start_time_ma_dur = as.numeric(as.Date("2020-05-15") - set_start_date()), # Change in length of stay post May 15, 2020 ("Medical Advances")
                       ##### Masking (NOT USED) ####
                       mask_cr_mat = NULL, # estimated contact reduction matrix including masking
                       start_time_mask_policy = Inf, # the starting time of the masking policy 
                       ##### Timestep ####
                       timestep = 0.05) # note, time step must be 1 or less
{
  
  for (x in c(1:length(epi_states_index))) {
    assign(names(epi_states_index)[x], epi_states_index[[x]])
  }
  
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
  
  #### Contact reduction matrix under SAH
  sip_flat <- c(t(sip_cr_mat))
  sip_contact_reduction <- lapply(c(1:length(sip_flat)), function(x) {
    sip_flat[x] * hash_table$m[[x]]
  })
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
  
  parms$prob_icu_death_no_bed <- c(prob_icu_death_no_bed)
  parms$prop_inf_die <- c(prop_inf_die)
  parms$prop_asymptomatic <- c(prop_asymptomatic)
  
  ###########################################
  #### Converting and rescale parameters ####
  ###########################################
  ## This part includes code in convert_parms.R and 
  ## parameter calculation in covid_19_model_function_v3.R
  
  timestep <- parms$timestep
  
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
  
  return(parms)
}
