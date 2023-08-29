#' @title Whether a dial back measure is greater than the threshold
#' 
#' @param icu_demand ICU beds at step t
#' @param avg_case2x Average case doubling time at step t
#' @param avg_hosp_7days Average new hospitalization in the past 7 days
#' @param prevI Prevalence at a time step
#' @param strategy A list contains the parameters for the cycling or dialback strategy
#' 
#' @export
get_db_criteria <- function(icu_demand, avg_case2x, avg_hosp_7days, prevI, 
                            strategy) {
  dial_back_measures <- list(icu = icu_demand, 
                             case2x = if(avg_case2x > 0) { -avg_case2x } else { -Inf }, 
                             hosp = avg_hosp_7days, 
                             prev = prevI)
  trigger_thres <- strategy$trigger_thres
  trigger_thres[["case2x"]] <- -trigger_thres[["case2x"]]
  
  out_criteria <- (dial_back_measures[[strategy$trigger]] >= 
                     trigger_thres[[strategy$trigger]])
  
  return(out_criteria)
}

#' @title Updating the start time of the trigger strategy
#' 
#' @param tnext Next time step
#' @param criteria Criteria calculated beforehand
#' @param parms List of parameters
#' 
#' @export
update_start_time_trigger_strategy <- function(tnext, criteria, parms) {
  tmp_time <- tnext * criteria
  start_time_sip <- tmp_time * (tmp_time > 0) + parms$start_time_sip * (tmp_time == 0)
  
  if (start_time_sip > parms$start_time_sip) {
    grace_period <- parms$trigger_strategy$grace_period
    parms$start_time_sip <- start_time_sip
    parms$trigger_grace_start <- start_time_sip + grace_period
    parms$start_time_trigger <- start_time_sip
    parms$end_time_trigger <- start_time_sip + grace_period + parms$trigger_strategy$min_on_day
    parms$end_time_sip <- parms$end_time_trigger
    parms$end_time_behavior_change <- start_time_sip
    parms$sip_days_past_peak <- -Inf 
    parms$trigger_strategy$start_time <- c(parms$trigger_strategy$start_time, start_time_sip)
    
    original_reduction <- parms$behavior_change_contact_reduction
    contact_red <- parms$trigger_strategy$on_contact_red * 
      (1 - (parms$hash_table$m_60p + parms$hash_table$m_60p_w_young)) + 
      parms$trigger_strategy$on_contact_red_60p * (parms$hash_table$m_60p + parms$hash_table$m_60p_w_young)
    sip_diff <- (contact_red - original_reduction) / grace_period
    parms$sip_contact_reduction <- original_reduction 
    parms$sip_diff <- sip_diff

  }
  
  return(parms)
}

#' @title Updating the end time of the trigger strategy
#' 
#' @param tnext Next time step
#' @param parms List of parameters
#' 
#' @export
update_end_time_trigger_strategy <- function(tnext, 
                                             v_days_since_peak, 
                                             str_peak_type, 
                                             end_day_peak_trigger, 
                                             parms) {
  # Set end day 
  # Both the dialback meausures are lower than the criteria, 
  # and the NPI is changed because we passed the peak
  tmp_cri <- (v_days_since_peak[[str_peak_type]] > end_day_peak_trigger) & 
    (tnext >= parms$end_time_trigger)
  
  if (isTRUE(tmp_cri)) {
    parms$end_time_sip <- tnext + parms$trigger_strategy$grace_period
    parms$end_time_trigger <- NULL 
    parms$trigger_grace_end <- parms$end_time_sip
    parms$start_time_behavior_change <- parms$end_time_sip
    parms$end_time_behavior_change <- Inf
    
    parms$trigger_strategy$end_time <- c(parms$trigger_strategy$end_time, parms$end_time_sip)
    
  } else {
    parms$end_time_sip <- tnext + 1
  }
  
  return(parms)
}
