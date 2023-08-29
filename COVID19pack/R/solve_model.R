#' @title Solve Model
#' 
#' takes initial conditions, sequence of times, model function, 
#' and parameter list and calculates the state of the model at each day timestep
#' 
#' @param times vector of times to calculate model state at
#' @param func function which generates differences between model states, \code{\link{covid_19_main_function}}
#' @param parms list of parameters
#' @param return_full default = 0 for integer day output; set to 1 if returning full output matrix by time step
#' 
#' @return 
#' matrix of raw results, each row is 1 timestep, each column is a different model state
#' 
#' @export
solve_model <- function(times, func, parms, return_full = 0, return_incidence = FALSE) {
  
  # parameters
  timestep <- parms$timestep
  v_strat_status <- parms$v_strat_status
  v_exp_str <- parms$v_exp_str
  v_asym_inf_str <- parms$v_asym_inf_str
  v_inf_str <- parms$v_inf_str
  y <- parms$init_vec
  tot_pop <- parms$N
  str_peak_type <- parms$str_peak_type
  max_time <- max(times)
  parms$max_time <- max_time
  
  index <- parms$ls_ix$index
  ie_str <- parms$ls_ix$ie_str
  
  # traversing through the time sequence
  model_state <- y
  y["Time"] <- 1.0
  
  # initialize the output matrix
  name_v_strat_status <- names(v_strat_status)
  
  out <- matrix(0, nrow = ((length(times) - 1) * timestep + 1), ncol = length(names(y)) + length(v_strat_status))
  colnames(out) <- c(names(y), name_v_strat_status)
  
  out_full <- matrix(0, nrow = (length(times)), ncol = length(names(y)) + length(v_strat_status))
  colnames(out_full) <- c(names(y), name_v_strat_status)
  
  # add the first row to the matrix
  out[1, ] <- out_full[1, ] <- c(y, 0 * v_strat_status)
  
  # End of existing interventions
  end_time_existing_intervention <- as.numeric(as.Date("2020-09-01") - set_start_date())
  
  # vector of days since peak
  v_days_since_peak <- c("infections" = -Inf, "hospitalizations" = -Inf, "deaths" = -Inf, 
                         "new_cases" = -Inf, "new_sym_cases" = -Inf)
  
  # track dial back strategy if parms$dialback_strategy != NULL
  if (!is.null(parms$trigger_strategy)) {
    end_day_peak_trigger <- parms$trigger_strategy$days_past_peak
    parms$trigger_strategy$start_time <- c()
    parms$trigger_strategy$end_time <- c()
  }
  
  # dial-back measures 
  v_new_cases <- rep(0, length(times))
  v_new_cases_daily <- rep(0, max_time)
  v_new_sym_cases <- rep(0, length(times))
  v_new_sym_cases_daily <- rep(0, max_time)
  v_new_hosp <- rep(0, length(times))
  v_new_hosp_daily <- rep(0, max_time) 
  v_case2x_time <- rep(0, max_time) # case doubling time
  
  # indices for calculating infections, hospitalizations, and deaths
  ind_inf <- parms$hash_table$v_ind_inf
  ind_hosp <- parms$hash_table$v_ind_hosp
  ind_icu <- parms$hash_table$icu_index
  ind_death <- parms$hash_table$v_ind_death
  
  if (isTRUE(return_incidence)) { # creating a list of incidence records
    ls_inc <- list()
    
    tmp <- names(y)[parms$hash_table$v_AI_ind1]
    ls_inc$new_AI <- matrix(0, nrow = (length(times)), ncol = length(tmp))
    colnames(ls_inc$new_AI) <- tmp
    
    tmp <- names(y)[parms$hash_table$v_I_ind1]
    ls_inc$new_I <- matrix(0, nrow = (length(times)), ncol = length(tmp))
    colnames(ls_inc$new_I) <- tmp
    
    tmp <- names(y)[parms$hash_table$v_QAI_ind1]
    ls_inc$new_QAI <- matrix(0, nrow = (length(times)), ncol = length(tmp))
    colnames(ls_inc$new_QAI) <- tmp
    
    tmp <- names(y)[parms$hash_table$v_QI_ind1]
    ls_inc$new_QI <- matrix(0, nrow = (length(times)), ncol = length(tmp))
    colnames(ls_inc$new_QI) <- tmp
    
    tmp <- names(y)[parms$hash_table$v_H_ind]
    ls_inc$new_H <- matrix(0, nrow = (length(times)), ncol = length(tmp))
    colnames(ls_inc$new_H) <- tmp
    
    tmp <- names(y)[parms$hash_table$v_ICU_ind]
    ls_inc$new_ICU <- matrix(0, nrow = (length(times)), ncol = length(tmp))
    colnames(ls_inc$new_ICU) <- tmp
  }
  
  for (i in 1:(length(times) - 1)) {
    t <- times[i]
    
    # for each step, call the model function to get the differentials
    v_strat_active_tlast <- out_full[i,  name_v_strat_status]
    
    model_state <- as.vector(model_state)
    ls_res <- func(t, model_state, parms, v_days_since_peak, v_strat_active_tlast)
    parms$mixing_matrix <- ls_res$mixing_matrix
    new_state <- ls_res$d
    v_strat_active <- ls_res$v_strat_active
    new_state <- as.vector(new_state)
    new_cases <- ls_res$new_cases
    new_sym_cases <- ls_res$new_sym_cases
    v_new_cases[i] <- new_cases
    v_new_sym_cases[i] <- new_sym_cases
    new_hosp <- ls_res$new_hosp
    v_new_hosp[i] <- new_hosp
    
    if (isTRUE(return_incidence)) {
      for (x in c(1:length(ls_inc))) {
        tmp_name <- names(ls_inc)[x]
        ls_inc[[tmp_name]][i, ] <- ls_res$ls_incidence[[tmp_name]]        
      }
    }
    
    if (parms$reduction_modified == 1) {
      parms$reduction_modified <- 0
    }
    
    # update the model state and add time 
    model_state <- model_state + new_state
    output_row <- model_state
    output_row["Time"] <- times[i + 1]
    
    out_full[i + 1, ] <- c(output_row, v_strat_active)
    
    # update abbreviated output if at an integer time point
    tnext <- times[i + 1]
    if (tnext == as.integer(tnext)) {
      # add state to the output matrix
      out[tnext, ] <- c(output_row, v_strat_active)
      v_new_cases_daily[tnext] <- sum(v_new_cases[(times > (tnext - 1)) & (times <= tnext)])
      v_new_sym_cases_daily[tnext] <- sum(v_new_sym_cases[(times > (tnext - 1)) & (times <= tnext)])
      v_new_hosp_daily[tnext] <- sum(v_new_hosp[(times > (tnext - 1)) & (times <= tnext)])
      db_time <- (tnext - (tnext - 1)) / log2(v_new_cases_daily[tnext] / v_new_cases_daily[tnext - 1])
      v_case2x_time[tnext] <- db_time
      
      # Update counter of days since peak
      # differences in infections, hospitalizations, and deaths
      if (tnext == 1) {
        v_days_since_peak["infections"] <- 0
        v_days_since_peak["hospitalizations"] <- 0
        daily_deaths <- diff(out[1, ind_death])
        v_days_since_peak["deaths"] <- 0
        v_days_since_peak["new_cases"] <- 0
        v_days_since_peak["new_sym_cases"] <- 0
      } else  {
        
        v_tmp <- rep(1, tnext)
        
        if (!is.null(parms$trigger_strategy) & (tnext >= end_time_existing_intervention)) {
          tstart <- max(c(parms$start_time_sip,  
                          end_time_existing_intervention))
          v_tmp <- rep(0, tnext)
          v_tmp[tstart:tnext] <- 1
        }
        
        tmp_row_sum_inf <- rowSums(out[1:tnext, ind_inf] * v_tmp)
        tmp_row_sum_hosp <- rowSums(out[1:tnext, ind_hosp] * v_tmp)
        daily_deaths <- diff(rowSums(out[1:tnext, ind_death]))
        tmp_new_cases <- v_new_cases_daily[1:tnext] * v_tmp
        tmp_new_sym_cases <- v_new_sym_cases_daily[1:tnext] * v_tmp
        
        v_days_since_peak["infections"] <- tnext - which.max(tmp_row_sum_inf)
        v_days_since_peak["hospitalizations"] <- tnext - which.max(tmp_row_sum_hosp)
        v_days_since_peak["deaths"] <- tnext - which.max(c(0, daily_deaths) * v_tmp)
        v_days_since_peak["new_cases"] <- tnext - which.max(tmp_new_cases)
        v_days_since_peak["new_sym_cases"] <- tnext - which.max(tmp_new_sym_cases)
        
        # Dial back measures
        if (tnext >= 7) {
          tnext_minus7 <- tnext - 7 # this is to go back to the records 7 days prior
          icu_demand <- sum(out[tnext_minus7, ind_icu])
          tmp_ix <- (tnext_minus7 - 7):tnext_minus7
          tmp_ix <- tmp_ix[tmp_ix > 0]
          avg_case2x <- mean(v_case2x_time[tmp_ix])
          tmp_ix <- (tnext_minus7 - 7):tnext_minus7
          tmp_ix <- tmp_ix[tmp_ix > 0]
          avg_hosp_7days <- sum(v_new_hosp_daily[tmp_ix]) / tot_pop # 7-day total new admission in hospitals
          prevI <- sum(out[tnext_minus7, ind_inf]) / tot_pop
        }
      }
      
      ##### Trigger strategy: Begins #####
      ## Change to defined age-independent contact reduction for trigger period
      if (!is.null(parms$trigger_strategy) & (tnext >= end_time_existing_intervention)) {
          parms$reduction_modified <- 1
          parms$behavior_change_contact_reduction <- matrix(parms$trigger_strategy$off_contact_red,
                                                            ncol = 9, nrow = 9)
          parms$behavior_change_contact_reduction[7:9,7:9] <- parms$trigger_strategy$off_contact_red_60p
      }
      
      if (!is.null(parms$trigger_strategy)) {
        if (tnext >= end_time_existing_intervention) {
          
          # using the dial back measure
          tmp_cri <- get_db_criteria(icu_demand, avg_case2x, avg_hosp_7days, prevI, 
                                     parms$trigger_strategy)
          
          # Set starting day
          if (is.null(parms$start_time_trigger) & isTRUE(tmp_cri)) { 
            parms <- update_start_time_trigger_strategy(tnext, tmp_cri, parms)
            # print(paste0("trigger on:"," day ", tnext))
          }
          
          if (!is.null(parms$trigger_grace_start)) {
            if (tnext < parms$trigger_grace_start) {
              parms$sip_contact_reduction <- parms$sip_contact_reduction + parms$sip_diff
              parms$reduction_modified <- 1
            }
          }
          
          # Set end day
          if (tnext == parms$end_time_sip & !is.null(parms$end_time_trigger)) {
            # either passed the peak or lower than the dial-back measure
            parms <- update_end_time_trigger_strategy(tnext, 
                                                      v_days_since_peak,
                                                      str_peak_type,
                                                      end_day_peak_trigger, 
                                                      parms)
          }
          
          if (!is.null(parms$trigger_grace_end)) {
            if (tnext < parms$trigger_grace_end) {
              parms$sip_contact_reduction <- parms$sip_contact_reduction - parms$sip_diff
              parms$reduction_modified <- 1
            }
            if (tnext == parms$trigger_grace_end) {
              parms$start_time_trigger <- NULL
              parms$trigger_grace_end <- NULL
            }
            
          }
        }
      }
      ##### Trigger strategy: Ends ##### 
      
      ##### Masking: Begins #####
      if (tnext == parms$start_time_mask_policy) {
        parms$behavior_change_contact_reduction <- parms$mask_cr_mat
        parms$behavior_change_contact_reduction[7:9, ] <- 0
        parms$behavior_change_contact_reduction[, 7:9] <- 0
        
        parms$sixty_plus_contact_reduction <- parms$mask_cr_mat
        parms$sixty_plus_contact_reduction[1:6, 1:6] <- 0
        
        parms$reduction_modified <- 1
      }
      
      ##### Masking: Ends #####
      
      
      ##### update prop_hosp
      if (any(tnext == parms$prop_hosp_change_day)) {
        tmp_ix <- which(tnext == parms$prop_hosp_change_day)
        parms$prop_hosp <- parms$prop_hosp_ls[[tmp_ix]]
      }
      
      ##### update prop_ICU 
      if (any(tnext == parms$prop_ICU_change_day)) {
        tmp_ix <- which(tnext == parms$prop_ICU_change_day)
        parms$prop_ICU <- parms$prop_ICU_ls[[tmp_ix]]
      }
      
    }
  }
  
  # return the matrix
  out[, "Time"] <- out[, "Time"] - 1
  out_full[, "Time"] <- out_full[, "Time"] - 1
  
  ls_out <- list() 
  if(return_full == 1) {
    ls_out$out_full <- out_full
  } else {
    ls_out$out <- out
  }
  
  if (isTRUE(return_incidence)) {
    ls_inc <- lapply(ls_inc, function(x) {
      x <- cbind(x, "Time" = times)
    })
    ls_out$ls_incidence <- ls_inc
  }
  
  if (!is.null(parms$trigger_strategy)) {
    ls_out$trigger_strategy <- parms$trigger_strategy
  }
  
  return(ls_out)
}
