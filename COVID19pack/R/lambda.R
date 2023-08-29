#' @title Modify mixing matrix via contact reduction 
#' @param t time step.
#' @param parms a list of parameters.
#' @param v_days_since_peak a vector of days since the peak of the measure specified in `str_peak_type`. 
#' @v_strat_active_tlast a vector of active NPIs at each time step. This vector is used to determine whether there are any changes in NPI that requires an updated calculation of the mixing matrix. 
#' @export
modify_mixing_matrix <- function(t, parms, v_days_since_peak, v_strat_active_tlast) {
  stsd <- parms$start_time_social_distancing
  stsc <- parms$start_time_school_closure
  stsip <- parms$start_time_sip
  st60p <- parms$start_time_60plus_distancing
  stbec <- parms$start_time_behavior_change
  et60p <- parms$end_time_60plus_distancing
  etsd <- parms$end_time_social_distancing
  etsc <- parms$end_time_school_closure
  etsip <- parms$end_time_sip
  etbec <- parms$end_time_behavior_change
  et_peak_sd <- parms$social_distancing_days_past_peak
  et_peak_sip <- parms$sip_days_past_peak
  et_peak_60p <- parms$sixty_plus_days_past_peak
  et_peak_sc <- parms$school_closure_days_past_peak
  et_peak_bec <- parms$behavior_change_days_past_peak
  str_peak_type <- parms$str_peak_type
  time <- t
  v_strat_status <- parms$v_strat_status
  mixing_matrix <- parms$mixing_matrix
  m_60p <- parms$hash_table[["m_60p"]]
  m_60p_w_young <- parms$hash_table[["m_60p_w_young"]]
  m_school <- parms$hash_table[["m_school"]]
  m_10to30 <- parms$hash_table[["m_10to30"]]
  reduction_modified <- parms$reduction_modified
  
  #### Creating operation factors 
  ### These are all dummy variables 
  social_distancing_op <- ((time >= stsd) &
                             (time < etsd | et_peak_sd > v_days_since_peak[str_peak_type]) & 
                             (v_strat_status["sd"] == 1)) * 1
  
  sip_op <- ((time >= stsip) & 
               (time < etsip | et_peak_sip > v_days_since_peak[str_peak_type]) & 
               (v_strat_status["sip"] == 1)) * 1
  
  school_closure_op <- ((time >= stsc) & 
                          (time < etsc | et_peak_sc > v_days_since_peak[str_peak_type]) & 
                          (v_strat_status["sc"] == 1)) * 1
  
  social_distancing_60_op <- ((time >= st60p) & 
                                (time < et60p | time < (st60p + et_peak_60p) | v_days_since_peak[str_peak_type] < et_peak_60p) & 
                                (v_strat_status["sd60p"] == 1)) * 1
  
  behavior_change_op <- ((time >= stbec) &
                           (time < etbec | time < (stbec + et_peak_bec) | v_days_since_peak[str_peak_type] < et_peak_bec) && 
                           (v_strat_status["bec"] == 1)) * 1
  
  v_op <- c("sd" = social_distancing_op[[1]], "sip" = sip_op[[1]], "sc" = school_closure_op[[1]], 
            "sd60p" = social_distancing_60_op[[1]], "bec" = behavior_change_op[[1]])
  if ((sip_op + social_distancing_op) == 2) {
    social_distancing_op <- 0
  }
  
  if (any(v_op != v_strat_active_tlast) | (reduction_modified == 1)) {
    mixing_matrix_o <- parms$mixing_matrix_o
    
    mixing_matrix <- sip_op * mixing_matrix_o * (1 - parms$sip_contact_reduction) + # if sip_op = 1, this mixing_matrix calculation would be activated
      social_distancing_op * mixing_matrix_o * (1 - parms$social_distancing_contact_reduction) + # if social_distancing_op = 1, this mixing_matrix calculation would be activated
      (1 - (sip_op + social_distancing_op)) * mixing_matrix_o * # if both sip_op and social_distancing_op are 0, this part of mixing_matrix calculation would be activated
      (1 - (school_closure_op * parms$school_closure_contact_reduction + # if both sip_op and social_distancing_op are 0 but school_closure_op = 1 and/or social_distancing_60_op = 1
              social_distancing_60_op * parms$sixty_plus_contact_reduction)) * 
      ((1 - behavior_change_op) +  
         behavior_change_op * ((school_closure_op * m_school + 
                                  social_distancing_60_op * (m_60p + m_60p_w_young)) +
      (1 - (school_closure_op * m_school + social_distancing_60_op * (m_60p + m_60p_w_young))) * 
        (1 - parms$behavior_change_contact_reduction)))
    # print("=====================")
    # print(t)
    # print(v_op)
    # print(mixing_matrix)
  }
  
  return(list(mixing_matrix = mixing_matrix, 
              v_op = v_op))
}

#' @title Calculate Lambda
#'
#' takes contact matrix, current state of model, and parameter set and calculates the number of contacts
#' each age group have with infected people
#'
#' @param mixing_matrix matrix which has the number of contacts age group i has with age group j e.g. mixing_matrix[1,2]
#' is the number of contacts age group 1 has with age group 2, 
#' @param init_pop current state of the model
#' @param parms a list of parameters
#'
#' @return
#' a vector which has the number of contacts a person in each age group has with an infected person each day
#'
#' @export
calculate_lambda <- function(mixing_matrix, init_pop, parms) {

  quarantine_contact_reduction <- parms$quarantine_contact_reduction
  N_by_age <- parms$N_by_age
  nag <- parms$n_age_groups
  ncg <- parms$n_co_groups
  
  v_ind_inf_by_age <- parms$hash_table[["v_ind_inf_by_age"]]
  v_ind_quar_inf_by_age <- parms$hash_table[["v_ind_quar_inf_by_age"]]
  v_ind_D_by_age <- parms$hash_table[["v_ind_D_by_age"]]
  v_ind_quar_noninf_by_age <- parms$hash_table[["v_ind_quar_noninf_by_age"]]
  
  # Calculate number of people in each age group who are infected and living
  # Number of infected by age
  n_inf_by_age <- colSums(matrix(init_pop[v_ind_inf_by_age], ncol = nag)) + 
    (1 - quarantine_contact_reduction) * colSums(matrix(init_pop[v_ind_quar_inf_by_age], ncol = nag))
  
  # Number of alive by age: remove non-infectious quarantined people who won't be circulating
  n_alive_by_age <- N_by_age - 
    colSums(matrix(init_pop[v_ind_D_by_age], ncol = nag)) - 
    quarantine_contact_reduction * colSums(matrix(init_pop[v_ind_quar_noninf_by_age], ncol = nag))

  # Calculate prevalence of infectious disease in each age group
  p_inf_by_age <- n_inf_by_age / n_alive_by_age

  # Calculate how many contacts a person in each age group has with other age groups and how many are infectious
  n_inf_contacts <- apply(mixing_matrix, 1, function(x) sum(x * p_inf_by_age))
  # Expand the contact vector by the number of comorbidity groups for each age group
  n_inf_contacts <- rep(n_inf_contacts, each = ncg) 

  return(lambda = n_inf_contacts)
}
