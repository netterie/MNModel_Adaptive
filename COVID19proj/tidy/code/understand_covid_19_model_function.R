# ****************************************************************************
# Code to understand covid_19_model_function(), located in COVID19pack/R 
# and called from solve_model.R (its alias is "func"). 
# ****************************************************************************

# ----------------------------------------------------------------------------
# Set up to step inside covid_19_model_function
# ----------------------------------------------------------------------------
library(tidyverse)
this_folder <- file.path("COVID19proj", "tidy", "code")

# Set up parms
source(file.path("COVID19proj", "tests", "test_scenarios.R"))
parms <- test_scenarios(nsamp=1,runmodel=FALSE)

# Initialize inputs to covid_19_model_function (alias: func) based on solve_model.R
# The call from solve_model is:
# ls_res <- func(t, model_state, parms, v_days_since_peak, v_strat_active_tlast)

# ----------------------------------------------------------------------------
# t
# ----------------------------------------------------------------------------
# time step 1
t           <- 1
# ----------------------------------------------------------------------------
# v_model_state
# ----------------------------------------------------------------------------
# from solve_model.R:
# model_state is the initial population vector, whereas y is exactly the same save
# for a 481st entry tracking the time step
y             <- parms$init_vec
model_state   <- y
model_state   <- as.vector(model_state)
y["Time"]     <- 1.0

v_model_state <- model_state

# ----------------------------------------------------------------------------
# parms
# ----------------------------------------------------------------------------
# as-is

# ----------------------------------------------------------------------------
# v_days_since_peak
# ----------------------------------------------------------------------------
v_days_since_peak <- c("infections" = -Inf, "hospitalizations" = -Inf, "deaths" = -Inf, 
                       "new_cases" = -Inf, "new_sym_cases" = -Inf)


# ----------------------------------------------------------------------------
# v_strat_active_tlast
# ----------------------------------------------------------------------------
# as in Scenarios.R
times <- seq(1, 365, by = parms$timestep)

# This is v_strat_status:
#    sd   sip    sc sd60p   bec 
#     1     1     1     1     1 
v_strat_status <- parms$v_strat_status
timestep       <- parms$timestep
# initialize the output matrices
# "out" has 1 row per day (365 total)
# (not time step) and one column per compartment PLUS one column for Time PLUS 
# 5 more columns, one for each entry in v_strat_status
# "out_full" is similar but has one row per time step
name_v_strat_status <- names(v_strat_status)
out <- matrix(0, nrow = ((length(times) - 1) * timestep + 1), ncol = length(names(y)) + length(v_strat_status))
colnames(out) <- c(names(y), name_v_strat_status)
out_full <- matrix(0, nrow = (length(times)), ncol = length(names(y)) + length(v_strat_status))
colnames(out_full) <- c(names(y), name_v_strat_status)
# add the first row to the matrix: init_vec, Time=1, and all v_strat_status entries = 0
out[1, ] <- out_full[1, ] <- c(y, 0 * v_strat_status)
# This next line gets updated inside the time loop: return the 
# most recent state of the v_strat_status entries
v_strat_active_tlast <- out_full[t,  name_v_strat_status]
#v_strat_active_tlast
#sd   sip    sc sd60p   bec 
#0     0     0     0     0 

# Call to function:
# ls_res <- func(t, model_state, parms, v_days_since_peak, v_strat_active_tlast)

# ----------------------------------------------------------------------------
#' @title Model Function for MN covid-19 model
#'
#' function contains a list of discrete differential equations
#' uses given parameter and timestep to calculate difference in each timestep at specified time
#'
#' states are
#' S-sucsceptible
#' E-pre-infectious, E1-Em
#' AI- asymptomatic infectious, AI1-AIn
#' I-infectious, I1-In
#' H- person in hospital
#' ICU-person in ICU
#' R-recovered
#' D-deaths
#'
#' @param t current time
#' @param v_model_state current state of model, i.e. number of people in each state, comorbidity, and age group
#' @param parms list of parameters
#' @param v_days_since_peak vector of the number of days since the last peak in infections, hospitalizations, and deaths
#'
#' @return vector of the difference in each model state at specified time step
#'
#' @export
# covid_19_model_function <- function(t, v_model_state, parms, v_days_since_peak, v_strat_active_tlast) {


# ----------------------------------------------------------------------------
# Initialize all parameters
# ----------------------------------------------------------------------------
  beta <- parms$beta
  beta_after <- parms$beta_after
  N <- parms$N
  nag <- parms$n_age_groups
  neg <- parms$n_epi_groups
  ncg <- parms$n_co_groups
  nes <- parms$n_exposed_states
  nis <- parms$n_infected_states
  str_peak_type <- parms$str_peak_type
  prop_hosp <- parms$prop_hosp
  hosp_mul <- parms$hosp_mul
  prop_ICU <- parms$prop_ICU
  time <- t
  prop_asym <- parms$prop_asymptomatic
  prop_inf_die <- parms$prop_inf_die
  prob_icu_death <- parms$prob_icu_death[["null"]]
  prob_hosp_death <-parms$prob_hosp_death[["null"]]
  prob_icu_death_no_bed <- parms$prob_icu_death_no_bed
  # Testing and quarantine states (deactivated)
  start_time_testing <- Inf # testing/quarantine is turned off
  p_spec_test <- parms$p_spec_test
  p_sens_test <- parms$p_sens_test
  actual_frac_SEAR_tested <- 0
  actual_frac_I_tested <- 0
  n_tests_per_day <- 0
  quarantine_contact_reduction <- parms$quarantine_contact_reduction
  p_I_seek_test <- parms$p_I_seek_test
  p_nonI_seek_test <- parms$p_nonI_seek_test
  p_trans_exp <- parms$p_trans_exp
  p_trans_inf <- parms$p_trans_inf
  p_quar_exit <- parms$p_quar_exit
  n_tests_per_day_hosp <- 0
  n_tests_per_day_icu <- 0
  # Hospialization and ICU use
  p_icu_exit_nobed <- parms$p_icu_exit_nobed
  p_h_exit <- parms$p_h_exit[["null"]]
  p_icu_exit <- parms$p_icu_exit[["null"]]
  hash_table <- parms$hash_table

# ----------------------------------------------------------------------------
# Update time-varying parameters that change based on calendar time
# ----------------------------------------------------------------------------

  # Set beta after September 1, 2020
  if (t>=163){
    beta <- beta_after
    prop_hosp <- prop_hosp * hosp_mul
  }

  ## Check if Medical Advances have occurred
  if(t >= parms$start_time_ma_mort){
    prob_hosp_death <- parms$prob_hosp_death[["MedAdv"]]
    prob_icu_death <- parms$prob_icu_death[["MedAdv"]]
  }

  if(t >= parms$start_time_ma_dur){
    p_h_exit <- parms$p_h_exit[["MedAdv"]]
    p_icu_exit <- parms$p_icu_exit[["MedAdv"]]
  }

# ----------------------------------------------------------------------------
# HEALTHCARE UTILIZATION: TESTING
# ----------------------------------------------------------------------------
# This section didn't get used, so I am not including it here



# ----------------------------------------------------------------------------
# TRANSMISSION: FORCE OF INFECTION and exiting S due to infection
# (Could have an external force of infection, e.g. travel or animal contact:
#  add entry rate to E/I from external force)
# ----------------------------------------------------------------------------

  ## Calculate number of contacts per day by age group based on current social distancing practices

  ## using current time, check if social distancing practices are in place
  ls_mixing_op <- modify_mixing_matrix(t, parms, v_days_since_peak, v_strat_active_tlast)
  # mixing_matrix is a n_agegroup x n_agegroup matrix (not "flattened" or anything)
  mixing_matrix <- ls_mixing_op$mixing_matrix

  # Call calculate_lambda() to get force of infection for current time step
  # lambda is a vector of length 18, one entry for each age-comorbidity group,
  # which is the length of v_S_ind.
  # calculate_lambda creates a vector of length 9, one entry per age group,
  # and then copies each entry using rep(...,each=ncg) 
  lambda <- calculate_lambda(mixing_matrix, v_model_state, parms)
  
    # Step inside calculate_lambda: 
    init_pop=v_model_state

  ##initilize vector of differences##
  d <- rep(0, length(parms$init_vec))

  ## dS
  v_S_ind <- hash_table$v_S_ind
  v_QS_ind <- hash_table$v_QS_ind

  # beta: (scalar = 0.02526144) probability that contact with an infectious person will cause an infection
  # lambda: (by compartment) # of infectious contacts
  # actual_frac_SEAR_tested: (scalar = 0)
  # p_quar_exit (scalar = 0.003565059)
  # p_spec_test (scalar = 0.99)
  d[v_S_ind] <- -beta * lambda *
    (1 - actual_frac_SEAR_tested * (1 - p_spec_test)) * v_model_state[v_S_ind] + # exiting to E1
    p_quar_exit * v_model_state[v_QS_ind] - # incoming
    actual_frac_SEAR_tested * (1 - p_spec_test) * v_model_state[v_S_ind] # exiting to QS

# ----------------------------------------------------------------------------
# DISEASE PROGRESSION: E-->AI or SI
# How long until you get symptoms (mild/moderate/severe) and (maybe) test
# positive
# ----------------------------------------------------------------------------

  ## dE1...dEm
  v_E_ind <- hash_table$v_E_ind
  # E1
  v_E_ind1 <- hash_table$v_E_ind1
  d[v_E_ind1] <- beta * lambda *
    (1 - actual_frac_SEAR_tested * (1 - p_spec_test)) * v_model_state[v_S_ind] - # incoming
    (1 - actual_frac_SEAR_tested * (1 - p_spec_test)) * p_trans_exp * v_model_state[v_E_ind1] - # outflow to E2 (if there are multiple exposed states) or to AI1 or I1 (if no other exposed states)
    actual_frac_SEAR_tested * (1 - p_spec_test) * v_model_state[v_E_ind1] # outflow to QE1

  # E2,...Em
  if (nes > 1) {
    v_E_ind_rm_last <- hash_table$v_E_ind_rm_last
    v_E_ind_rm_first <- hash_table$v_E_ind_rm_first
    d[v_E_ind_rm_first] <-
      p_trans_exp * (1 - actual_frac_SEAR_tested * (1 - p_spec_test)) * v_model_state[v_E_ind_rm_last] - # incoming E_(i-1) --> E_i
      (1 - actual_frac_SEAR_tested * (1 - p_spec_test)) * p_trans_exp * v_model_state[v_E_ind_rm_first] - # exiting E_i --> E_(i+1)
      actual_frac_SEAR_tested * (1 - p_spec_test) * v_model_state[v_E_ind_rm_first] # outflow to the corresponding QE_i state
  }

  ## dAI1...dAIn
  v_AI_ind <- hash_table$v_AI_ind
  # AI1
  v_AI_ind1 <- hash_table$v_AI_ind1
  v_E_ind_last <- hash_table$v_E_ind_last
  #x <- prop_asym[cag + 1, agg + 1] * p_trans_exp * v_model_state[v_E_ind[nes]] - p_trans_inf * v_model_state[v_AI_ind[1]]
  d[v_AI_ind1] <- prop_asym * (1 - actual_frac_SEAR_tested * (1 - p_spec_test)) *
    p_trans_exp * v_model_state[v_E_ind_last] - # inflow from the last E state
    p_trans_inf * (1 - p_sens_test * actual_frac_SEAR_tested) * v_model_state[v_AI_ind1] - # outflow to AI_2
    p_sens_test * actual_frac_SEAR_tested * v_model_state[v_AI_ind1] # outflow to QAI1

  new_AI <- prop_asym * (1 - actual_frac_SEAR_tested * (1 - p_spec_test)) *
    p_trans_exp * v_model_state[v_E_ind_last]

  # AI2,...AIn
  if (nis > 1) {
    v_AI_ind_rm_last <- hash_table$v_AI_ind_rm_last
    v_AI_ind_rm_first <- hash_table$v_AI_ind_rm_first
    d[v_AI_ind_rm_first] <-
      p_trans_inf * (1 - p_sens_test * actual_frac_SEAR_tested) * (v_model_state[v_AI_ind_rm_last] - # incoming I_(i-1) --> I_i
         v_model_state[v_AI_ind_rm_first]) - # exiting I_i --> I_(i+1)
      p_sens_test * actual_frac_SEAR_tested * v_model_state[v_AI_ind_rm_first] # outflow to
  }

  ## dI1...dIn
  v_I_ind <- hash_table$v_I_ind
  # I1
  v_I_ind1 <- hash_table$v_I_ind1
  d[v_I_ind1] <- (1 - prop_asym) * (1 - actual_frac_SEAR_tested * (1 - p_spec_test)) *
    p_trans_exp * v_model_state[v_E_ind_last] - # inflow from the last E_n state
    p_trans_inf * (1 - p_sens_test * actual_frac_I_tested) * v_model_state[v_I_ind1] - # outflow to I2
    p_sens_test * actual_frac_I_tested * v_model_state[v_I_ind1] # outflow to QI1

  new_I <- (1 - prop_asym) * (1 - actual_frac_SEAR_tested * (1 - p_spec_test)) *
    p_trans_exp * v_model_state[v_E_ind_last]

  # I2,...In
  if (nis > 1) {
    v_I_ind_rm_last <- hash_table$v_I_ind_rm_last
    v_I_ind_rm_first <- hash_table$v_I_ind_rm_first
    d[v_I_ind_rm_first] <-
      p_trans_inf * (1 - p_sens_test * actual_frac_I_tested) *
      (v_model_state[v_I_ind_rm_last] - # incoming I_(i-1) --> I_i
         v_model_state[v_I_ind_rm_first]) - # exiting I_i --> I_(i+1)
      p_sens_test * actual_frac_I_tested * v_model_state[v_I_ind_rm_first]
  }

  ## dQS
  v_QS_ind <- hash_table$v_QS_ind
  d[v_QS_ind] <- actual_frac_SEAR_tested * (1 - p_spec_test) * v_model_state[v_S_ind] - # inflow from S
    p_quar_exit * v_model_state[v_QS_ind] - # outflow to S
    beta * (1 - quarantine_contact_reduction) * lambda * v_model_state[v_QS_ind] # outflow to QE

  new_case <- sum(actual_frac_SEAR_tested * (1 - p_spec_test) * v_model_state[v_S_ind])

  ## dQE1...dQEm
  v_QE_ind <- hash_table$v_QE_ind
  # QE1
  v_QE_ind1 <- hash_table$v_QE_ind1
  d[v_QE_ind1] <- beta * (1 - quarantine_contact_reduction) * lambda * v_model_state[v_QS_ind] + #inflow from QS
    actual_frac_SEAR_tested * (1 - p_spec_test) * v_model_state[v_E_ind1] - # inflow from E1
    p_trans_exp * v_model_state[v_QE_ind1] # outflow to QE2

  new_case <- new_case + sum(actual_frac_SEAR_tested * (1 - p_spec_test) * v_model_state[v_E_ind1])

  # QE2,...QEm
  if (nes > 1) {
    v_QE_ind_rm_last <- hash_table$v_QE_ind_rm_last
    v_QE_ind_rm_first <- hash_table$v_QE_ind_rm_first
    d[v_QE_ind_rm_first] <-
      p_trans_exp * v_model_state[v_QE_ind_rm_last] - # incoming E_(i-1) --> E_i
      p_trans_exp * v_model_state[v_QE_ind_rm_first] + # exiting E_i --> E_(i+1)
      actual_frac_SEAR_tested * (1 - p_spec_test) * v_model_state[v_E_ind_rm_first]

    new_case <- new_case + sum(actual_frac_SEAR_tested * (1 - p_spec_test) * v_model_state[v_E_ind_rm_first])
  }

  ## dQAI1...dQAIn
  v_QAI_ind <- hash_table$v_QAI_ind
  # QAI1
  v_QAI_ind1 <- hash_table$v_QAI_ind1
  v_QE_ind_last <- hash_table$v_QE_ind_last
  d[v_QAI_ind1] <- prop_asym * p_trans_exp * v_model_state[v_QE_ind_last] - # inflow from QE_n
    p_trans_inf * v_model_state[v_QAI_ind1] + #outflow to next QAI_i
    p_sens_test * actual_frac_SEAR_tested * v_model_state[v_AI_ind1] #inflow from AI1

  new_case <- new_case + sum(p_sens_test * actual_frac_SEAR_tested * v_model_state[v_AI_ind1])
  new_QAI <- prop_asym * p_trans_exp * v_model_state[v_QE_ind_last]

  # QAI2,...QAIn
  if (nis > 1) {
    v_QAI_ind_rm_last <- hash_table$v_QAI_ind_rm_last
    v_QAI_ind_rm_first <- hash_table$v_QAI_ind_rm_first

    d[v_QAI_ind_rm_first] <-
      p_trans_inf * (v_model_state[v_QAI_ind_rm_last] - # incoming I_(i-1) --> I_i (inflow from previous AI_i state)
                       v_model_state[v_QAI_ind_rm_first]) + # exiting I_i --> I_(i+1) (outflow to next AI state)
      p_sens_test * actual_frac_SEAR_tested * v_model_state[v_AI_ind_rm_first]

    new_case <- new_case + sum(p_sens_test * actual_frac_SEAR_tested * v_model_state[v_AI_ind_rm_first])
  }

  ## dQI1...dQIn
  # v_QI_ind <- v_ind$v_QI_ind
  v_QI_ind <- hash_table$v_QI_ind
  # QI1
  v_QI_ind1 <- hash_table$v_QI_ind1
  d[v_QI_ind1] <- (1 - prop_asym) * p_trans_exp * v_model_state[v_QE_ind_last] - # inflow from QE1
    p_trans_inf * v_model_state[v_QI_ind1] + # outflow to next QI_i
    p_sens_test * actual_frac_I_tested * v_model_state[v_I_ind1] # inflow from I1

  new_case <- new_case + sum(p_sens_test * actual_frac_I_tested * v_model_state[v_I_ind1])
  new_QI <- (1 - prop_asym) * p_trans_exp * v_model_state[v_QE_ind_last]

  # QI2,...QIn
  if (nis > 1) {
    v_QI_ind_rm_last <- hash_table$v_QI_ind_rm_last
    v_QI_ind_rm_first <- hash_table$v_QI_ind_rm_first

    d[v_QI_ind_rm_first] <-
      p_trans_inf * (v_model_state[v_QI_ind_rm_last] - # incoming I_(i-1) --> I_i
                       v_model_state[v_QI_ind_rm_first]) + # exiting I_i --> I_(i+1)
      p_sens_test * actual_frac_I_tested * v_model_state[v_I_ind_rm_first]

    new_case <- new_case + sum(p_sens_test * actual_frac_I_tested * v_model_state[v_I_ind_rm_first])
  }

# Add exiting to H or ICU here

# ----------------------------------------------------------------------------
# HEALTHCARE UTILIZATION: I--> HOSPITAL OR ICU (once you get there)
# For a different disease, there could be different disease severities that
# exit to H/ICU at different rates
# ----------------------------------------------------------------------------

  ## H (hospitalization)
  v_H_ind <- hash_table$v_H_ind
  v_I_ind_last <- hash_table$v_I_ind_last
  v_QI_ind_last <- hash_table$v_QI_ind_last
  d[v_H_ind] <-
    p_trans_inf * prop_hosp * (1 - prop_ICU) *
    (v_model_state[v_I_ind_last] + v_model_state[v_QI_ind_last]) - # incoming from infected state
    p_h_exit * v_model_state[v_H_ind]

  new_H <- p_trans_inf * prop_hosp * (1 - prop_ICU) *
    (v_model_state[v_I_ind_last] + v_model_state[v_QI_ind_last])


  ## Check for ICU over-capacity

  # Calculate total number of people in ICU
  v_ind_icu <- hash_table[["icu_index"]]
  n_icu_t <- sum(v_model_state[v_ind_icu])
  if (is.na(n_icu_t)) {
    stop(paste0("ICU numbers are negative!!: time= ", time))
  }

  # Calculate the proportion of people who need an ICU but cannot access one
  p_icu_overflow <- 0
  if (n_icu_t > parms$n_icu_beds) {
    p_icu_overflow <- (n_icu_t - parms$n_icu_beds) / n_icu_t
  }

  ## ICU
  v_ICU_ind <- hash_table$v_ICU_ind
  # weighted average of exit probabilities (due to potential ICU overwhelm)
  p_icu_exit_overall <- (1 - p_icu_overflow) * p_icu_exit +
    p_icu_overflow * p_icu_exit_nobed
  # dICU
  d[v_ICU_ind] <-
    p_trans_inf * prop_hosp * prop_ICU *
    (v_model_state[v_I_ind_last] + v_model_state[v_QI_ind_last]) - # incoming from infected state
    p_icu_exit_overall * v_model_state[v_ICU_ind] # ICU overall probability of exiting

  new_case <- new_case + sum(p_trans_inf * prop_hosp * v_model_state[v_I_ind_last])
  new_ICU <- p_trans_inf * prop_hosp * prop_ICU *
    (v_model_state[v_I_ind_last] + v_model_state[v_QI_ind_last])

  ## CH (cumulative hospitalizations)
  v_CH_ind <- hash_table$v_CH_ind
  # dCH
  d[v_CH_ind] <-
    p_trans_inf * (v_model_state[v_I_ind_last] + v_model_state[v_QI_ind_last]) *
    prop_hosp * (1 - prop_ICU)

# ----------------------------------------------------------------------------
# Combine R, D or R-with-complications into a "Post-Infection Module"
# Can they go back into S? Depends on prophylaxis, vaccination...
# ----------------------------------------------------------------------------

# ----------------------------------------------------------------------------
# RECOVERY: Weave into AI, SI, H and ICU exit sections?
# ----------------------------------------------------------------------------

  ## R (recovery)
  v_R_ind <- hash_table$v_R_ind
  # weighted average of recovery flow (due to potential ICU overwhelm)
  p_icu2rec_overall <- (1 - p_icu_overflow) * p_icu_exit * (1 - prob_icu_death) +
    p_icu_overflow * p_icu_exit_nobed * (1 - prob_icu_death_no_bed)

  # dR
  v_AI_ind_last <- hash_table$v_AI_ind_last
  d[v_R_ind] <- (1 - prop_inf_die) * (1 - prop_hosp) * p_trans_inf * v_model_state[v_I_ind_last] + # incoming directly from infectious, no hospitalization
    p_trans_inf * (v_model_state[v_AI_ind_last])  # incoming from aysmptomatic infection

  # dRD
  v_RD_ind <- hash_table$v_RD_ind
  v_QAI_ind_last <- hash_table$v_QAI_ind_last
  d[v_RD_ind] <- p_icu2rec_overall * v_model_state[v_ICU_ind] + # recovering from ICU
    (1 - prob_hosp_death) * p_h_exit * v_model_state[v_H_ind] + #recovering from hospital
    (1 - prop_inf_die) * p_trans_inf * (v_model_state[v_QI_ind_last]) * (1 - prop_hosp) + #recovering from QI
    p_trans_inf * v_model_state[v_QAI_ind_last] #recovering from qai

# ----------------------------------------------------------------------------
# DEATH: Weave into SI, H and ICU exit sections?
# ----------------------------------------------------------------------------

  ## D (death)
  v_D_ind <- hash_table$v_D_ind
  # dD
  d[v_D_ind] <-
    p_icu_exit * (1 - p_icu_overflow) * prob_icu_death * v_model_state[v_ICU_ind] + # incoming from ICU, with bed
    p_icu_exit_nobed * p_icu_overflow * prob_icu_death_no_bed * v_model_state[v_ICU_ind] + # incoming from ICU, without bed
    prop_inf_die * p_trans_inf * (v_model_state[v_I_ind_last] + v_model_state[v_QI_ind_last]) * (1 - prop_hosp) + # incoming from Infectious state
    prob_hosp_death * p_h_exit * v_model_state[v_H_ind]

  ## HD (home death)
  v_HD_ind <- hash_table$v_HD_ind
  # dHD
  d[v_HD_ind] <-
    prop_inf_die * p_trans_inf * (v_model_state[v_I_ind_last] + v_model_state[v_QI_ind_last]) *
    (1 - prop_hosp)

# ----------------------------------------------------------------------------
# COLLECT EPI STATS
# ----------------------------------------------------------------------------


  ls_incidence <- list(new_AI = new_AI, new_I = new_I,
                       new_QAI = new_QAI, new_QI = new_QI,
                       new_H = new_H, new_ICU = new_ICU)

  return(list("d" = d,
              "new_cases" = new_case,
              "new_sym_cases" = sum(new_I) + sum(new_QI),
              "new_hosp" = sum(new_H) + sum(new_ICU),
              "v_strat_active" = ls_mixing_op$v_op,
              "mixing_matrix" = mixing_matrix,
              "ls_incidence" = ls_incidence))
}
#' @title Model Function for MN covid-19 model
#'
#' function contains a list of discrete differential equations
#' uses given parameter and timestep to calculate difference in each timestep at specified time
#'
#' states are
#' S-sucsceptible
#' E-pre-infectious, E1-Em
#' AI- asymptomatic infectious, AI1-AIn
#' I-infectious, I1-In
#' H- person in hospital
#' ICU-person in ICU
#' R-recovered
#' D-deaths
#'
#' @param t current time
#' @param v_model_state current state of model, i.e. number of people in each state, comorbidity, and age group
#' @param parms list of parameters
#' @param v_days_since_peak vector of the number of days since the last peak in infections, hospitalizations, and deaths
#'
#' @return vector of the difference in each model state at specified time step
#'
#' @export
covid_19_model_function <- function(t, v_model_state, parms, v_days_since_peak, v_strat_active_tlast) {

# ----------------------------------------------------------------------------
# Initialize all parameters
# ----------------------------------------------------------------------------
  beta <- parms$beta
  beta_after <- parms$beta_after
  N <- parms$N
  nag <- parms$n_age_groups
  neg <- parms$n_epi_groups
  ncg <- parms$n_co_groups
  nes <- parms$n_exposed_states
  nis <- parms$n_infected_states
  str_peak_type <- parms$str_peak_type
  prop_hosp <- parms$prop_hosp
  hosp_mul <- parms$hosp_mul
  prop_ICU <- parms$prop_ICU
  time <- t
  prop_asym <- parms$prop_asymptomatic
  prop_inf_die <- parms$prop_inf_die
  prob_icu_death <- parms$prob_icu_death[["null"]]
  prob_hosp_death <-parms$prob_hosp_death[["null"]]
  prob_icu_death_no_bed <- parms$prob_icu_death_no_bed
  # Testing and quarantine states (deactivated)
  start_time_testing <- Inf # testing/quarantine is turned off
  p_spec_test <- parms$p_spec_test
  p_sens_test <- parms$p_sens_test
  actual_frac_SEAR_tested <- 0
  actual_frac_I_tested <- 0
  n_tests_per_day <- 0
  quarantine_contact_reduction <- parms$quarantine_contact_reduction
  p_I_seek_test <- parms$p_I_seek_test
  p_nonI_seek_test <- parms$p_nonI_seek_test
  p_trans_exp <- parms$p_trans_exp
  p_trans_inf <- parms$p_trans_inf
  p_quar_exit <- parms$p_quar_exit
  n_tests_per_day_hosp <- 0
  n_tests_per_day_icu <- 0
  # Hospialization and ICU use
  p_icu_exit_nobed <- parms$p_icu_exit_nobed
  p_h_exit <- parms$p_h_exit[["null"]]
  p_icu_exit <- parms$p_icu_exit[["null"]]
  hash_table <- parms$hash_table

# ----------------------------------------------------------------------------
# Update time-varying parameters that change based on calendar time
# ----------------------------------------------------------------------------

  # Set beta after September 1, 2020
  if (t>=163){
    beta <- beta_after
    prop_hosp <- prop_hosp * hosp_mul
  }

  ## Check if Medical Advances have occurred
  if(t >= parms$start_time_ma_mort){
    prob_hosp_death <- parms$prob_hosp_death[["MedAdv"]]
    prob_icu_death <- parms$prob_icu_death[["MedAdv"]]
  }

  if(t >= parms$start_time_ma_dur){
    p_h_exit <- parms$p_h_exit[["MedAdv"]]
    p_icu_exit <- parms$p_icu_exit[["MedAdv"]]
  }

# ----------------------------------------------------------------------------
# HEALTHCARE UTILIZATION: TESTING
# ----------------------------------------------------------------------------

# Update time-varying parameters that change based on model state
  ## Increase in Testing- Fill in actual date
  if(t > start_time_testing) {
    v_ind_s <- hash_table[["s_index"]]
    v_ind_exp <- hash_table[["exp_index"]]
    v_ind_inf <- hash_table[["inf_index"]]
    v_ind_asym <- hash_table[["asym_index"]]
    v_ind_rec <- hash_table[["rec_index"]]
    v_ind_icu <- hash_table[["icu_index"]]
    v_ind_hosp <- hash_table[["hosp_index"]]

    n_s_pop <- sum(v_model_state[v_ind_s])
    n_exp_pop <- sum(v_model_state[v_ind_exp])
    n_inf_pop <- sum(v_model_state[v_ind_inf])
    n_asym_pop <- sum(v_model_state[v_ind_asym])
    n_rec_pop <- sum(v_model_state[v_ind_rec])
    n_icu_pop <- sum(v_model_state[v_ind_icu])
    n_hosp_pop <- sum(v_model_state[v_ind_hosp])
    n_nonI_pop <- n_s_pop + n_exp_pop + n_asym_pop + n_rec_pop
    # adjust number of tests per day based on number needed in ICU and Hospital

    n_tests_per_day <- max(c(0, n_tests_per_day -
                               n_tests_per_day_hosp * n_hosp_pop -
                               n_tests_per_day_icu * n_icu_pop))

    # calculate total demand for testing
    n_test_demand <- p_I_seek_test * n_inf_pop + p_nonI_seek_test * n_nonI_pop

    # calculate actual proportions tested, applying testing capacity constraint
    if (n_test_demand > n_tests_per_day) {
      # probability of being tested in I state, given proportion of testing demand that can be filled
      p_I_tested <- p_I_seek_test * (n_tests_per_day/n_test_demand)
      # probability of being tested in non-I state, given proportion of testing demand that can be filled
      p_nonI_tested <- p_nonI_seek_test * (n_tests_per_day/n_test_demand)

    } else {
      # if demand does not exceed capacity, all tests are filled
      p_I_tested <- p_I_seek_test
      p_nonI_tested <- p_nonI_seek_test
    }
    actual_frac_SEAR_tested <- p_nonI_tested
    actual_frac_I_tested <- p_I_tested
  }



# ----------------------------------------------------------------------------
# TRANSMISSION: FORCE OF INFECTION and exiting S due to infection
# (Could have an external force of infection, e.g. travel or animal contact:
#  add entry rate to E/I from external force)
# ----------------------------------------------------------------------------

  ## Calculate number of contacts per day by age group based on current social distancing practices

  ## using current time, check if social distancing practices are in place
  ls_mixing_op <- modify_mixing_matrix(t, parms, v_days_since_peak, v_strat_active_tlast)
  mixing_matrix <- ls_mixing_op$mixing_matrix

  # Call calculate_lambda() to get force of infection for current time step
  lambda <- calculate_lambda(mixing_matrix, v_model_state, parms)

  ##initilize vector of differences##
  d <- rep(0, length(parms$init_vec))

  ## dS
  v_S_ind <- hash_table$v_S_ind
  v_QS_ind <- hash_table$v_QS_ind

  d[v_S_ind] <- -beta * lambda *
    (1 - actual_frac_SEAR_tested * (1 - p_spec_test)) * v_model_state[v_S_ind] + # exiting to E1
    p_quar_exit * v_model_state[v_QS_ind] - # incoming
    actual_frac_SEAR_tested * (1 - p_spec_test) * v_model_state[v_S_ind] # exiting to QS

# ----------------------------------------------------------------------------
# DISEASE PROGRESSION: E-->AI or SI
# How long until you get symptoms (mild/moderate/severe) and (maybe) test
# positive
# ----------------------------------------------------------------------------

  ## dE1...dEm
  v_E_ind <- hash_table$v_E_ind
  # E1
  v_E_ind1 <- hash_table$v_E_ind1
  d[v_E_ind1] <- beta * lambda *
    (1 - actual_frac_SEAR_tested * (1 - p_spec_test)) * v_model_state[v_S_ind] - # incoming
    (1 - actual_frac_SEAR_tested * (1 - p_spec_test)) * p_trans_exp * v_model_state[v_E_ind1] - # outflow to E2 (if there are multiple exposed states) or to AI1 or I1 (if no other exposed states)
    actual_frac_SEAR_tested * (1 - p_spec_test) * v_model_state[v_E_ind1] # outflow to QE1

  # E2,...Em
  if (nes > 1) {
    v_E_ind_rm_last <- hash_table$v_E_ind_rm_last
    v_E_ind_rm_first <- hash_table$v_E_ind_rm_first
    d[v_E_ind_rm_first] <-
      p_trans_exp * (1 - actual_frac_SEAR_tested * (1 - p_spec_test)) * v_model_state[v_E_ind_rm_last] - # incoming E_(i-1) --> E_i
      (1 - actual_frac_SEAR_tested * (1 - p_spec_test)) * p_trans_exp * v_model_state[v_E_ind_rm_first] - # exiting E_i --> E_(i+1)
      actual_frac_SEAR_tested * (1 - p_spec_test) * v_model_state[v_E_ind_rm_first] # outflow to the corresponding QE_i state
  }

  ## dAI1...dAIn
  v_AI_ind <- hash_table$v_AI_ind
  # AI1
  v_AI_ind1 <- hash_table$v_AI_ind1
  v_E_ind_last <- hash_table$v_E_ind_last
  #x <- prop_asym[cag + 1, agg + 1] * p_trans_exp * v_model_state[v_E_ind[nes]] - p_trans_inf * v_model_state[v_AI_ind[1]]
  d[v_AI_ind1] <- prop_asym * (1 - actual_frac_SEAR_tested * (1 - p_spec_test)) *
    p_trans_exp * v_model_state[v_E_ind_last] - # inflow from the last E state
    p_trans_inf * (1 - p_sens_test * actual_frac_SEAR_tested) * v_model_state[v_AI_ind1] - # outflow to AI_2
    p_sens_test * actual_frac_SEAR_tested * v_model_state[v_AI_ind1] # outflow to QAI1

  new_AI <- prop_asym * (1 - actual_frac_SEAR_tested * (1 - p_spec_test)) *
    p_trans_exp * v_model_state[v_E_ind_last]

  # AI2,...AIn
  if (nis > 1) {
    v_AI_ind_rm_last <- hash_table$v_AI_ind_rm_last
    v_AI_ind_rm_first <- hash_table$v_AI_ind_rm_first
    d[v_AI_ind_rm_first] <-
      p_trans_inf * (1 - p_sens_test * actual_frac_SEAR_tested) * (v_model_state[v_AI_ind_rm_last] - # incoming I_(i-1) --> I_i
         v_model_state[v_AI_ind_rm_first]) - # exiting I_i --> I_(i+1)
      p_sens_test * actual_frac_SEAR_tested * v_model_state[v_AI_ind_rm_first] # outflow to
  }

  ## dI1...dIn
  v_I_ind <- hash_table$v_I_ind
  # I1
  v_I_ind1 <- hash_table$v_I_ind1
  d[v_I_ind1] <- (1 - prop_asym) * (1 - actual_frac_SEAR_tested * (1 - p_spec_test)) *
    p_trans_exp * v_model_state[v_E_ind_last] - # inflow from the last E_n state
    p_trans_inf * (1 - p_sens_test * actual_frac_I_tested) * v_model_state[v_I_ind1] - # outflow to I2
    p_sens_test * actual_frac_I_tested * v_model_state[v_I_ind1] # outflow to QI1

  new_I <- (1 - prop_asym) * (1 - actual_frac_SEAR_tested * (1 - p_spec_test)) *
    p_trans_exp * v_model_state[v_E_ind_last]

  # I2,...In
  if (nis > 1) {
    v_I_ind_rm_last <- hash_table$v_I_ind_rm_last
    v_I_ind_rm_first <- hash_table$v_I_ind_rm_first
    d[v_I_ind_rm_first] <-
      p_trans_inf * (1 - p_sens_test * actual_frac_I_tested) *
      (v_model_state[v_I_ind_rm_last] - # incoming I_(i-1) --> I_i
         v_model_state[v_I_ind_rm_first]) - # exiting I_i --> I_(i+1)
      p_sens_test * actual_frac_I_tested * v_model_state[v_I_ind_rm_first]
  }

  ## dQS
  v_QS_ind <- hash_table$v_QS_ind
  d[v_QS_ind] <- actual_frac_SEAR_tested * (1 - p_spec_test) * v_model_state[v_S_ind] - # inflow from S
    p_quar_exit * v_model_state[v_QS_ind] - # outflow to S
    beta * (1 - quarantine_contact_reduction) * lambda * v_model_state[v_QS_ind] # outflow to QE

  new_case <- sum(actual_frac_SEAR_tested * (1 - p_spec_test) * v_model_state[v_S_ind])

  ## dQE1...dQEm
  v_QE_ind <- hash_table$v_QE_ind
  # QE1
  v_QE_ind1 <- hash_table$v_QE_ind1
  d[v_QE_ind1] <- beta * (1 - quarantine_contact_reduction) * lambda * v_model_state[v_QS_ind] + #inflow from QS
    actual_frac_SEAR_tested * (1 - p_spec_test) * v_model_state[v_E_ind1] - # inflow from E1
    p_trans_exp * v_model_state[v_QE_ind1] # outflow to QE2

  new_case <- new_case + sum(actual_frac_SEAR_tested * (1 - p_spec_test) * v_model_state[v_E_ind1])

  # QE2,...QEm
  if (nes > 1) {
    v_QE_ind_rm_last <- hash_table$v_QE_ind_rm_last
    v_QE_ind_rm_first <- hash_table$v_QE_ind_rm_first
    d[v_QE_ind_rm_first] <-
      p_trans_exp * v_model_state[v_QE_ind_rm_last] - # incoming E_(i-1) --> E_i
      p_trans_exp * v_model_state[v_QE_ind_rm_first] + # exiting E_i --> E_(i+1)
      actual_frac_SEAR_tested * (1 - p_spec_test) * v_model_state[v_E_ind_rm_first]

    new_case <- new_case + sum(actual_frac_SEAR_tested * (1 - p_spec_test) * v_model_state[v_E_ind_rm_first])
  }

  ## dQAI1...dQAIn
  v_QAI_ind <- hash_table$v_QAI_ind
  # QAI1
  v_QAI_ind1 <- hash_table$v_QAI_ind1
  v_QE_ind_last <- hash_table$v_QE_ind_last
  d[v_QAI_ind1] <- prop_asym * p_trans_exp * v_model_state[v_QE_ind_last] - # inflow from QE_n
    p_trans_inf * v_model_state[v_QAI_ind1] + #outflow to next QAI_i
    p_sens_test * actual_frac_SEAR_tested * v_model_state[v_AI_ind1] #inflow from AI1

  new_case <- new_case + sum(p_sens_test * actual_frac_SEAR_tested * v_model_state[v_AI_ind1])
  new_QAI <- prop_asym * p_trans_exp * v_model_state[v_QE_ind_last]

  # QAI2,...QAIn
  if (nis > 1) {
    v_QAI_ind_rm_last <- hash_table$v_QAI_ind_rm_last
    v_QAI_ind_rm_first <- hash_table$v_QAI_ind_rm_first

    d[v_QAI_ind_rm_first] <-
      p_trans_inf * (v_model_state[v_QAI_ind_rm_last] - # incoming I_(i-1) --> I_i (inflow from previous AI_i state)
                       v_model_state[v_QAI_ind_rm_first]) + # exiting I_i --> I_(i+1) (outflow to next AI state)
      p_sens_test * actual_frac_SEAR_tested * v_model_state[v_AI_ind_rm_first]

    new_case <- new_case + sum(p_sens_test * actual_frac_SEAR_tested * v_model_state[v_AI_ind_rm_first])
  }

  ## dQI1...dQIn
  # v_QI_ind <- v_ind$v_QI_ind
  v_QI_ind <- hash_table$v_QI_ind
  # QI1
  v_QI_ind1 <- hash_table$v_QI_ind1
  d[v_QI_ind1] <- (1 - prop_asym) * p_trans_exp * v_model_state[v_QE_ind_last] - # inflow from QE1
    p_trans_inf * v_model_state[v_QI_ind1] + # outflow to next QI_i
    p_sens_test * actual_frac_I_tested * v_model_state[v_I_ind1] # inflow from I1

  new_case <- new_case + sum(p_sens_test * actual_frac_I_tested * v_model_state[v_I_ind1])
  new_QI <- (1 - prop_asym) * p_trans_exp * v_model_state[v_QE_ind_last]

  # QI2,...QIn
  if (nis > 1) {
    v_QI_ind_rm_last <- hash_table$v_QI_ind_rm_last
    v_QI_ind_rm_first <- hash_table$v_QI_ind_rm_first

    d[v_QI_ind_rm_first] <-
      p_trans_inf * (v_model_state[v_QI_ind_rm_last] - # incoming I_(i-1) --> I_i
                       v_model_state[v_QI_ind_rm_first]) + # exiting I_i --> I_(i+1)
      p_sens_test * actual_frac_I_tested * v_model_state[v_I_ind_rm_first]

    new_case <- new_case + sum(p_sens_test * actual_frac_I_tested * v_model_state[v_I_ind_rm_first])
  }

# Add exiting to H or ICU here

# ----------------------------------------------------------------------------
# HEALTHCARE UTILIZATION: I--> HOSPITAL OR ICU (once you get there)
# For a different disease, there could be different disease severities that
# exit to H/ICU at different rates
# ----------------------------------------------------------------------------

  ## H (hospitalization)
  v_H_ind <- hash_table$v_H_ind
  v_I_ind_last <- hash_table$v_I_ind_last
  v_QI_ind_last <- hash_table$v_QI_ind_last
  d[v_H_ind] <-
    p_trans_inf * prop_hosp * (1 - prop_ICU) *
    (v_model_state[v_I_ind_last] + v_model_state[v_QI_ind_last]) - # incoming from infected state
    p_h_exit * v_model_state[v_H_ind]

  new_H <- p_trans_inf * prop_hosp * (1 - prop_ICU) *
    (v_model_state[v_I_ind_last] + v_model_state[v_QI_ind_last])


  ## Check for ICU over-capacity

  # Calculate total number of people in ICU
  v_ind_icu <- hash_table[["icu_index"]]
  n_icu_t <- sum(v_model_state[v_ind_icu])
  if (is.na(n_icu_t)) {
    stop(paste0("ICU numbers are negative!!: time= ", time))
  }

  # Calculate the proportion of people who need an ICU but cannot access one
  p_icu_overflow <- 0
  if (n_icu_t > parms$n_icu_beds) {
    p_icu_overflow <- (n_icu_t - parms$n_icu_beds) / n_icu_t
  }

  ## ICU
  v_ICU_ind <- hash_table$v_ICU_ind
  # weighted average of exit probabilities (due to potential ICU overwhelm)
  p_icu_exit_overall <- (1 - p_icu_overflow) * p_icu_exit +
    p_icu_overflow * p_icu_exit_nobed
  # dICU
  d[v_ICU_ind] <-
    p_trans_inf * prop_hosp * prop_ICU *
    (v_model_state[v_I_ind_last] + v_model_state[v_QI_ind_last]) - # incoming from infected state
    p_icu_exit_overall * v_model_state[v_ICU_ind] # ICU overall probability of exiting

  new_case <- new_case + sum(p_trans_inf * prop_hosp * v_model_state[v_I_ind_last])
  new_ICU <- p_trans_inf * prop_hosp * prop_ICU *
    (v_model_state[v_I_ind_last] + v_model_state[v_QI_ind_last])

  ## CH (cumulative hospitalizations)
  v_CH_ind <- hash_table$v_CH_ind
  # dCH
  d[v_CH_ind] <-
    p_trans_inf * (v_model_state[v_I_ind_last] + v_model_state[v_QI_ind_last]) *
    prop_hosp * (1 - prop_ICU)

# ----------------------------------------------------------------------------
# Combine R, D or R-with-complications into a "Post-Infection Module"
# Can they go back into S? Depends on prophylaxis, vaccination...
# ----------------------------------------------------------------------------

# ----------------------------------------------------------------------------
# RECOVERY: Weave into AI, SI, H and ICU exit sections?
# ----------------------------------------------------------------------------

  ## R (recovery)
  v_R_ind <- hash_table$v_R_ind
  # weighted average of recovery flow (due to potential ICU overwhelm)
  p_icu2rec_overall <- (1 - p_icu_overflow) * p_icu_exit * (1 - prob_icu_death) +
    p_icu_overflow * p_icu_exit_nobed * (1 - prob_icu_death_no_bed)

  # dR
  v_AI_ind_last <- hash_table$v_AI_ind_last
  d[v_R_ind] <- (1 - prop_inf_die) * (1 - prop_hosp) * p_trans_inf * v_model_state[v_I_ind_last] + # incoming directly from infectious, no hospitalization
    p_trans_inf * (v_model_state[v_AI_ind_last])  # incoming from aysmptomatic infection

  # dRD
  v_RD_ind <- hash_table$v_RD_ind
  v_QAI_ind_last <- hash_table$v_QAI_ind_last
  d[v_RD_ind] <- p_icu2rec_overall * v_model_state[v_ICU_ind] + # recovering from ICU
    (1 - prob_hosp_death) * p_h_exit * v_model_state[v_H_ind] + #recovering from hospital
    (1 - prop_inf_die) * p_trans_inf * (v_model_state[v_QI_ind_last]) * (1 - prop_hosp) + #recovering from QI
    p_trans_inf * v_model_state[v_QAI_ind_last] #recovering from qai

# ----------------------------------------------------------------------------
# DEATH: Weave into SI, H and ICU exit sections?
# ----------------------------------------------------------------------------

  ## D (death)
  v_D_ind <- hash_table$v_D_ind
  # dD
  d[v_D_ind] <-
    p_icu_exit * (1 - p_icu_overflow) * prob_icu_death * v_model_state[v_ICU_ind] + # incoming from ICU, with bed
    p_icu_exit_nobed * p_icu_overflow * prob_icu_death_no_bed * v_model_state[v_ICU_ind] + # incoming from ICU, without bed
    prop_inf_die * p_trans_inf * (v_model_state[v_I_ind_last] + v_model_state[v_QI_ind_last]) * (1 - prop_hosp) + # incoming from Infectious state
    prob_hosp_death * p_h_exit * v_model_state[v_H_ind]

  ## HD (home death)
  v_HD_ind <- hash_table$v_HD_ind
  # dHD
  d[v_HD_ind] <-
    prop_inf_die * p_trans_inf * (v_model_state[v_I_ind_last] + v_model_state[v_QI_ind_last]) *
    (1 - prop_hosp)

# ----------------------------------------------------------------------------
# COLLECT EPI STATS
# ----------------------------------------------------------------------------


  ls_incidence <- list(new_AI = new_AI, new_I = new_I,
                       new_QAI = new_QAI, new_QI = new_QI,
                       new_H = new_H, new_ICU = new_ICU)

  return(list("d" = d,
              "new_cases" = new_case,
              "new_sym_cases" = sum(new_I) + sum(new_QI),
              "new_hosp" = sum(new_H) + sum(new_ICU),
              "v_strat_active" = ls_mixing_op$v_op,
              "mixing_matrix" = mixing_matrix,
              "ls_incidence" = ls_incidence))
}
#' @title Model Function for MN covid-19 model
#'
#' function contains a list of discrete differential equations
#' uses given parameter and timestep to calculate difference in each timestep at specified time
#'
#' states are
#' S-sucsceptible
#' E-pre-infectious, E1-Em
#' AI- asymptomatic infectious, AI1-AIn
#' I-infectious, I1-In
#' H- person in hospital
#' ICU-person in ICU
#' R-recovered
#' D-deaths
#'
#' @param t current time
#' @param v_model_state current state of model, i.e. number of people in each state, comorbidity, and age group
#' @param parms list of parameters
#' @param v_days_since_peak vector of the number of days since the last peak in infections, hospitalizations, and deaths
#'
#' @return vector of the difference in each model state at specified time step
#'
#' @export
covid_19_model_function <- function(t, v_model_state, parms, v_days_since_peak, v_strat_active_tlast) {

# ----------------------------------------------------------------------------
# Initialize all parameters
# ----------------------------------------------------------------------------
  beta <- parms$beta
  beta_after <- parms$beta_after
  N <- parms$N
  nag <- parms$n_age_groups
  neg <- parms$n_epi_groups
  ncg <- parms$n_co_groups
  nes <- parms$n_exposed_states
  nis <- parms$n_infected_states
  str_peak_type <- parms$str_peak_type
  prop_hosp <- parms$prop_hosp
  hosp_mul <- parms$hosp_mul
  prop_ICU <- parms$prop_ICU
  time <- t
  prop_asym <- parms$prop_asymptomatic
  prop_inf_die <- parms$prop_inf_die
  prob_icu_death <- parms$prob_icu_death[["null"]]
  prob_hosp_death <-parms$prob_hosp_death[["null"]]
  prob_icu_death_no_bed <- parms$prob_icu_death_no_bed
  # Testing and quarantine states (deactivated)
  start_time_testing <- Inf # testing/quarantine is turned off
  p_spec_test <- parms$p_spec_test
  p_sens_test <- parms$p_sens_test
  actual_frac_SEAR_tested <- 0
  actual_frac_I_tested <- 0
  n_tests_per_day <- 0
  quarantine_contact_reduction <- parms$quarantine_contact_reduction
  p_I_seek_test <- parms$p_I_seek_test
  p_nonI_seek_test <- parms$p_nonI_seek_test
  p_trans_exp <- parms$p_trans_exp
  p_trans_inf <- parms$p_trans_inf
  p_quar_exit <- parms$p_quar_exit
  n_tests_per_day_hosp <- 0
  n_tests_per_day_icu <- 0
  # Hospialization and ICU use
  p_icu_exit_nobed <- parms$p_icu_exit_nobed
  p_h_exit <- parms$p_h_exit[["null"]]
  p_icu_exit <- parms$p_icu_exit[["null"]]
  hash_table <- parms$hash_table

# ----------------------------------------------------------------------------
# Update time-varying parameters that change based on calendar time
# ----------------------------------------------------------------------------

  # Set beta after September 1, 2020
  if (t>=163){
    beta <- beta_after
    prop_hosp <- prop_hosp * hosp_mul
  }

  ## Check if Medical Advances have occurred
  if(t >= parms$start_time_ma_mort){
    prob_hosp_death <- parms$prob_hosp_death[["MedAdv"]]
    prob_icu_death <- parms$prob_icu_death[["MedAdv"]]
  }

  if(t >= parms$start_time_ma_dur){
    p_h_exit <- parms$p_h_exit[["MedAdv"]]
    p_icu_exit <- parms$p_icu_exit[["MedAdv"]]
  }

# ----------------------------------------------------------------------------
# HEALTHCARE UTILIZATION: TESTING
# ----------------------------------------------------------------------------

# Update time-varying parameters that change based on model state
  ## Increase in Testing- Fill in actual date
  if(t > start_time_testing) {
    v_ind_s <- hash_table[["s_index"]]
    v_ind_exp <- hash_table[["exp_index"]]
    v_ind_inf <- hash_table[["inf_index"]]
    v_ind_asym <- hash_table[["asym_index"]]
    v_ind_rec <- hash_table[["rec_index"]]
    v_ind_icu <- hash_table[["icu_index"]]
    v_ind_hosp <- hash_table[["hosp_index"]]

    n_s_pop <- sum(v_model_state[v_ind_s])
    n_exp_pop <- sum(v_model_state[v_ind_exp])
    n_inf_pop <- sum(v_model_state[v_ind_inf])
    n_asym_pop <- sum(v_model_state[v_ind_asym])
    n_rec_pop <- sum(v_model_state[v_ind_rec])
    n_icu_pop <- sum(v_model_state[v_ind_icu])
    n_hosp_pop <- sum(v_model_state[v_ind_hosp])
    n_nonI_pop <- n_s_pop + n_exp_pop + n_asym_pop + n_rec_pop
    # adjust number of tests per day based on number needed in ICU and Hospital

    n_tests_per_day <- max(c(0, n_tests_per_day -
                               n_tests_per_day_hosp * n_hosp_pop -
                               n_tests_per_day_icu * n_icu_pop))

    # calculate total demand for testing
    n_test_demand <- p_I_seek_test * n_inf_pop + p_nonI_seek_test * n_nonI_pop

    # calculate actual proportions tested, applying testing capacity constraint
    if (n_test_demand > n_tests_per_day) {
      # probability of being tested in I state, given proportion of testing demand that can be filled
      p_I_tested <- p_I_seek_test * (n_tests_per_day/n_test_demand)
      # probability of being tested in non-I state, given proportion of testing demand that can be filled
      p_nonI_tested <- p_nonI_seek_test * (n_tests_per_day/n_test_demand)

    } else {
      # if demand does not exceed capacity, all tests are filled
      p_I_tested <- p_I_seek_test
      p_nonI_tested <- p_nonI_seek_test
    }
    actual_frac_SEAR_tested <- p_nonI_tested
    actual_frac_I_tested <- p_I_tested
  }



# ----------------------------------------------------------------------------
# TRANSMISSION: FORCE OF INFECTION and exiting S due to infection
# (Could have an external force of infection, e.g. travel or animal contact:
#  add entry rate to E/I from external force)
# ----------------------------------------------------------------------------

  ## Calculate number of contacts per day by age group based on current social distancing practices

  ## using current time, check if social distancing practices are in place
  ls_mixing_op <- modify_mixing_matrix(t, parms, v_days_since_peak, v_strat_active_tlast)
  mixing_matrix <- ls_mixing_op$mixing_matrix

  # Call calculate_lambda() to get force of infection for current time step
  lambda <- calculate_lambda(mixing_matrix, v_model_state, parms)

  ##initilize vector of differences##
  d <- rep(0, length(parms$init_vec))

  ## dS
  v_S_ind <- hash_table$v_S_ind
  v_QS_ind <- hash_table$v_QS_ind

  d[v_S_ind] <- -beta * lambda *
    (1 - actual_frac_SEAR_tested * (1 - p_spec_test)) * v_model_state[v_S_ind] + # exiting to E1
    p_quar_exit * v_model_state[v_QS_ind] - # incoming
    actual_frac_SEAR_tested * (1 - p_spec_test) * v_model_state[v_S_ind] # exiting to QS

# ----------------------------------------------------------------------------
# DISEASE PROGRESSION: E-->AI or SI
# How long until you get symptoms (mild/moderate/severe) and (maybe) test
# positive
# ----------------------------------------------------------------------------

  ## dE1...dEm
  v_E_ind <- hash_table$v_E_ind
  # E1
  v_E_ind1 <- hash_table$v_E_ind1
  d[v_E_ind1] <- beta * lambda *
    (1 - actual_frac_SEAR_tested * (1 - p_spec_test)) * v_model_state[v_S_ind] - # incoming
    (1 - actual_frac_SEAR_tested * (1 - p_spec_test)) * p_trans_exp * v_model_state[v_E_ind1] - # outflow to E2 (if there are multiple exposed states) or to AI1 or I1 (if no other exposed states)
    actual_frac_SEAR_tested * (1 - p_spec_test) * v_model_state[v_E_ind1] # outflow to QE1

  # E2,...Em
  if (nes > 1) {
    v_E_ind_rm_last <- hash_table$v_E_ind_rm_last
    v_E_ind_rm_first <- hash_table$v_E_ind_rm_first
    d[v_E_ind_rm_first] <-
      p_trans_exp * (1 - actual_frac_SEAR_tested * (1 - p_spec_test)) * v_model_state[v_E_ind_rm_last] - # incoming E_(i-1) --> E_i
      (1 - actual_frac_SEAR_tested * (1 - p_spec_test)) * p_trans_exp * v_model_state[v_E_ind_rm_first] - # exiting E_i --> E_(i+1)
      actual_frac_SEAR_tested * (1 - p_spec_test) * v_model_state[v_E_ind_rm_first] # outflow to the corresponding QE_i state
  }

  ## dAI1...dAIn
  v_AI_ind <- hash_table$v_AI_ind
  # AI1
  v_AI_ind1 <- hash_table$v_AI_ind1
  v_E_ind_last <- hash_table$v_E_ind_last
  #x <- prop_asym[cag + 1, agg + 1] * p_trans_exp * v_model_state[v_E_ind[nes]] - p_trans_inf * v_model_state[v_AI_ind[1]]
  d[v_AI_ind1] <- prop_asym * (1 - actual_frac_SEAR_tested * (1 - p_spec_test)) *
    p_trans_exp * v_model_state[v_E_ind_last] - # inflow from the last E state
    p_trans_inf * (1 - p_sens_test * actual_frac_SEAR_tested) * v_model_state[v_AI_ind1] - # outflow to AI_2
    p_sens_test * actual_frac_SEAR_tested * v_model_state[v_AI_ind1] # outflow to QAI1

  new_AI <- prop_asym * (1 - actual_frac_SEAR_tested * (1 - p_spec_test)) *
    p_trans_exp * v_model_state[v_E_ind_last]

  # AI2,...AIn
  if (nis > 1) {
    v_AI_ind_rm_last <- hash_table$v_AI_ind_rm_last
    v_AI_ind_rm_first <- hash_table$v_AI_ind_rm_first
    d[v_AI_ind_rm_first] <-
      p_trans_inf * (1 - p_sens_test * actual_frac_SEAR_tested) * (v_model_state[v_AI_ind_rm_last] - # incoming I_(i-1) --> I_i
         v_model_state[v_AI_ind_rm_first]) - # exiting I_i --> I_(i+1)
      p_sens_test * actual_frac_SEAR_tested * v_model_state[v_AI_ind_rm_first] # outflow to
  }

  ## dI1...dIn
  v_I_ind <- hash_table$v_I_ind
  # I1
  v_I_ind1 <- hash_table$v_I_ind1
  d[v_I_ind1] <- (1 - prop_asym) * (1 - actual_frac_SEAR_tested * (1 - p_spec_test)) *
    p_trans_exp * v_model_state[v_E_ind_last] - # inflow from the last E_n state
    p_trans_inf * (1 - p_sens_test * actual_frac_I_tested) * v_model_state[v_I_ind1] - # outflow to I2
    p_sens_test * actual_frac_I_tested * v_model_state[v_I_ind1] # outflow to QI1

  new_I <- (1 - prop_asym) * (1 - actual_frac_SEAR_tested * (1 - p_spec_test)) *
    p_trans_exp * v_model_state[v_E_ind_last]

  # I2,...In
  if (nis > 1) {
    v_I_ind_rm_last <- hash_table$v_I_ind_rm_last
    v_I_ind_rm_first <- hash_table$v_I_ind_rm_first
    d[v_I_ind_rm_first] <-
      p_trans_inf * (1 - p_sens_test * actual_frac_I_tested) *
      (v_model_state[v_I_ind_rm_last] - # incoming I_(i-1) --> I_i
         v_model_state[v_I_ind_rm_first]) - # exiting I_i --> I_(i+1)
      p_sens_test * actual_frac_I_tested * v_model_state[v_I_ind_rm_first]
  }

  ## dQS
  v_QS_ind <- hash_table$v_QS_ind
  d[v_QS_ind] <- actual_frac_SEAR_tested * (1 - p_spec_test) * v_model_state[v_S_ind] - # inflow from S
    p_quar_exit * v_model_state[v_QS_ind] - # outflow to S
    beta * (1 - quarantine_contact_reduction) * lambda * v_model_state[v_QS_ind] # outflow to QE

  new_case <- sum(actual_frac_SEAR_tested * (1 - p_spec_test) * v_model_state[v_S_ind])

  ## dQE1...dQEm
  v_QE_ind <- hash_table$v_QE_ind
  # QE1
  v_QE_ind1 <- hash_table$v_QE_ind1
  d[v_QE_ind1] <- beta * (1 - quarantine_contact_reduction) * lambda * v_model_state[v_QS_ind] + #inflow from QS
    actual_frac_SEAR_tested * (1 - p_spec_test) * v_model_state[v_E_ind1] - # inflow from E1
    p_trans_exp * v_model_state[v_QE_ind1] # outflow to QE2

  new_case <- new_case + sum(actual_frac_SEAR_tested * (1 - p_spec_test) * v_model_state[v_E_ind1])

  # QE2,...QEm
  if (nes > 1) {
    v_QE_ind_rm_last <- hash_table$v_QE_ind_rm_last
    v_QE_ind_rm_first <- hash_table$v_QE_ind_rm_first
    d[v_QE_ind_rm_first] <-
      p_trans_exp * v_model_state[v_QE_ind_rm_last] - # incoming E_(i-1) --> E_i
      p_trans_exp * v_model_state[v_QE_ind_rm_first] + # exiting E_i --> E_(i+1)
      actual_frac_SEAR_tested * (1 - p_spec_test) * v_model_state[v_E_ind_rm_first]

    new_case <- new_case + sum(actual_frac_SEAR_tested * (1 - p_spec_test) * v_model_state[v_E_ind_rm_first])
  }

  ## dQAI1...dQAIn
  v_QAI_ind <- hash_table$v_QAI_ind
  # QAI1
  v_QAI_ind1 <- hash_table$v_QAI_ind1
  v_QE_ind_last <- hash_table$v_QE_ind_last
  d[v_QAI_ind1] <- prop_asym * p_trans_exp * v_model_state[v_QE_ind_last] - # inflow from QE_n
    p_trans_inf * v_model_state[v_QAI_ind1] + #outflow to next QAI_i
    p_sens_test * actual_frac_SEAR_tested * v_model_state[v_AI_ind1] #inflow from AI1

  new_case <- new_case + sum(p_sens_test * actual_frac_SEAR_tested * v_model_state[v_AI_ind1])
  new_QAI <- prop_asym * p_trans_exp * v_model_state[v_QE_ind_last]

  # QAI2,...QAIn
  if (nis > 1) {
    v_QAI_ind_rm_last <- hash_table$v_QAI_ind_rm_last
    v_QAI_ind_rm_first <- hash_table$v_QAI_ind_rm_first

    d[v_QAI_ind_rm_first] <-
      p_trans_inf * (v_model_state[v_QAI_ind_rm_last] - # incoming I_(i-1) --> I_i (inflow from previous AI_i state)
                       v_model_state[v_QAI_ind_rm_first]) + # exiting I_i --> I_(i+1) (outflow to next AI state)
      p_sens_test * actual_frac_SEAR_tested * v_model_state[v_AI_ind_rm_first]

    new_case <- new_case + sum(p_sens_test * actual_frac_SEAR_tested * v_model_state[v_AI_ind_rm_first])
  }

  ## dQI1...dQIn
  # v_QI_ind <- v_ind$v_QI_ind
  v_QI_ind <- hash_table$v_QI_ind
  # QI1
  v_QI_ind1 <- hash_table$v_QI_ind1
  d[v_QI_ind1] <- (1 - prop_asym) * p_trans_exp * v_model_state[v_QE_ind_last] - # inflow from QE1
    p_trans_inf * v_model_state[v_QI_ind1] + # outflow to next QI_i
    p_sens_test * actual_frac_I_tested * v_model_state[v_I_ind1] # inflow from I1

  new_case <- new_case + sum(p_sens_test * actual_frac_I_tested * v_model_state[v_I_ind1])
  new_QI <- (1 - prop_asym) * p_trans_exp * v_model_state[v_QE_ind_last]

  # QI2,...QIn
  if (nis > 1) {
    v_QI_ind_rm_last <- hash_table$v_QI_ind_rm_last
    v_QI_ind_rm_first <- hash_table$v_QI_ind_rm_first

    d[v_QI_ind_rm_first] <-
      p_trans_inf * (v_model_state[v_QI_ind_rm_last] - # incoming I_(i-1) --> I_i
                       v_model_state[v_QI_ind_rm_first]) + # exiting I_i --> I_(i+1)
      p_sens_test * actual_frac_I_tested * v_model_state[v_I_ind_rm_first]

    new_case <- new_case + sum(p_sens_test * actual_frac_I_tested * v_model_state[v_I_ind_rm_first])
  }

# Add exiting to H or ICU here

# ----------------------------------------------------------------------------
# HEALTHCARE UTILIZATION: I--> HOSPITAL OR ICU (once you get there)
# For a different disease, there could be different disease severities that
# exit to H/ICU at different rates
# ----------------------------------------------------------------------------

  ## H (hospitalization)
  v_H_ind <- hash_table$v_H_ind
  v_I_ind_last <- hash_table$v_I_ind_last
  v_QI_ind_last <- hash_table$v_QI_ind_last
  d[v_H_ind] <-
    p_trans_inf * prop_hosp * (1 - prop_ICU) *
    (v_model_state[v_I_ind_last] + v_model_state[v_QI_ind_last]) - # incoming from infected state
    p_h_exit * v_model_state[v_H_ind]

  new_H <- p_trans_inf * prop_hosp * (1 - prop_ICU) *
    (v_model_state[v_I_ind_last] + v_model_state[v_QI_ind_last])


  ## Check for ICU over-capacity

  # Calculate total number of people in ICU
  v_ind_icu <- hash_table[["icu_index"]]
  n_icu_t <- sum(v_model_state[v_ind_icu])
  if (is.na(n_icu_t)) {
    stop(paste0("ICU numbers are negative!!: time= ", time))
  }

  # Calculate the proportion of people who need an ICU but cannot access one
  p_icu_overflow <- 0
  if (n_icu_t > parms$n_icu_beds) {
    p_icu_overflow <- (n_icu_t - parms$n_icu_beds) / n_icu_t
  }

  ## ICU
  v_ICU_ind <- hash_table$v_ICU_ind
  # weighted average of exit probabilities (due to potential ICU overwhelm)
  p_icu_exit_overall <- (1 - p_icu_overflow) * p_icu_exit +
    p_icu_overflow * p_icu_exit_nobed
  # dICU
  d[v_ICU_ind] <-
    p_trans_inf * prop_hosp * prop_ICU *
    (v_model_state[v_I_ind_last] + v_model_state[v_QI_ind_last]) - # incoming from infected state
    p_icu_exit_overall * v_model_state[v_ICU_ind] # ICU overall probability of exiting

  new_case <- new_case + sum(p_trans_inf * prop_hosp * v_model_state[v_I_ind_last])
  new_ICU <- p_trans_inf * prop_hosp * prop_ICU *
    (v_model_state[v_I_ind_last] + v_model_state[v_QI_ind_last])

  ## CH (cumulative hospitalizations)
  v_CH_ind <- hash_table$v_CH_ind
  # dCH
  d[v_CH_ind] <-
    p_trans_inf * (v_model_state[v_I_ind_last] + v_model_state[v_QI_ind_last]) *
    prop_hosp * (1 - prop_ICU)

# ----------------------------------------------------------------------------
# Combine R, D or R-with-complications into a "Post-Infection Module"
# Can they go back into S? Depends on prophylaxis, vaccination...
# ----------------------------------------------------------------------------

# ----------------------------------------------------------------------------
# RECOVERY: Weave into AI, SI, H and ICU exit sections?
# ----------------------------------------------------------------------------

  ## R (recovery)
  v_R_ind <- hash_table$v_R_ind
  # weighted average of recovery flow (due to potential ICU overwhelm)
  p_icu2rec_overall <- (1 - p_icu_overflow) * p_icu_exit * (1 - prob_icu_death) +
    p_icu_overflow * p_icu_exit_nobed * (1 - prob_icu_death_no_bed)

  # dR
  v_AI_ind_last <- hash_table$v_AI_ind_last
  d[v_R_ind] <- (1 - prop_inf_die) * (1 - prop_hosp) * p_trans_inf * v_model_state[v_I_ind_last] + # incoming directly from infectious, no hospitalization
    p_trans_inf * (v_model_state[v_AI_ind_last])  # incoming from aysmptomatic infection

  # dRD
  v_RD_ind <- hash_table$v_RD_ind
  v_QAI_ind_last <- hash_table$v_QAI_ind_last
  d[v_RD_ind] <- p_icu2rec_overall * v_model_state[v_ICU_ind] + # recovering from ICU
    (1 - prob_hosp_death) * p_h_exit * v_model_state[v_H_ind] + #recovering from hospital
    (1 - prop_inf_die) * p_trans_inf * (v_model_state[v_QI_ind_last]) * (1 - prop_hosp) + #recovering from QI
    p_trans_inf * v_model_state[v_QAI_ind_last] #recovering from qai

# ----------------------------------------------------------------------------
# DEATH: Weave into SI, H and ICU exit sections?
# ----------------------------------------------------------------------------

  ## D (death)
  v_D_ind <- hash_table$v_D_ind
  # dD
  d[v_D_ind] <-
    p_icu_exit * (1 - p_icu_overflow) * prob_icu_death * v_model_state[v_ICU_ind] + # incoming from ICU, with bed
    p_icu_exit_nobed * p_icu_overflow * prob_icu_death_no_bed * v_model_state[v_ICU_ind] + # incoming from ICU, without bed
    prop_inf_die * p_trans_inf * (v_model_state[v_I_ind_last] + v_model_state[v_QI_ind_last]) * (1 - prop_hosp) + # incoming from Infectious state
    prob_hosp_death * p_h_exit * v_model_state[v_H_ind]

  ## HD (home death)
  v_HD_ind <- hash_table$v_HD_ind
  # dHD
  d[v_HD_ind] <-
    prop_inf_die * p_trans_inf * (v_model_state[v_I_ind_last] + v_model_state[v_QI_ind_last]) *
    (1 - prop_hosp)

# ----------------------------------------------------------------------------
# COLLECT EPI STATS
# ----------------------------------------------------------------------------


  ls_incidence <- list(new_AI = new_AI, new_I = new_I,
                       new_QAI = new_QAI, new_QI = new_QI,
                       new_H = new_H, new_ICU = new_ICU)

  return(list("d" = d,
              "new_cases" = new_case,
              "new_sym_cases" = sum(new_I) + sum(new_QI),
              "new_hosp" = sum(new_H) + sum(new_ICU),
              "v_strat_active" = ls_mixing_op$v_op,
              "mixing_matrix" = mixing_matrix,
              "ls_incidence" = ls_incidence))
}