#' @title Process Output
#'
#' takes raw output from solve function and calculates cumulative infections, prevalent infections, daily deaths, cumulative deaths,
#' number of people in the hospital, number of people in the ICU for each timestep
#'
#' @param out matrix or list of raw output from solve function
#' @param parms list of parameters
#'
#' @return
#' the same matrix that was inputted but with columns at the start for each calculated output
#'
#' @import data.table
#' @export
process_output <- function(out, parms){
  
  neg <- parms$n_epi_groups
  nag <- parms$n_age_groups
  ncg <- parms$n_co_groups
  ts <- parms$timestep
  nes <- parms$n_exposed_states
  nis <- parms$n_infected_states
  
  index <- parms$ls_ix$index
  ie_str <- parms$ls_ix$ie_str
  ia <- parms$ls_ix$ia
  
  for(i in c(1:length(parms$epi_groups_ls))) {
    assign(names(parms$epi_groups_ls)[i], parms$epi_groups_ls[[i]])
  }

  chk <- !is.matrix(out) & (length(out) > 1)
  if (chk) {
    inc_data <- out$ls_incidence
  }
  out <- out$out
  
  ## extract prevalent infections by summing over E+AI+I+QE+QI+QAI+ICU+H 
  inf_ix <- index[ie_str %in% c(v_exp_str, v_asym_inf_str, v_inf_str,v_qe_str,v_qi_str,v_qai_str,
                                "H", "ICU")]
  prevalent_infections <- rowSums(out[, inf_ix])

  ## extract cumulative infections by summing over  E+AI+I+QE+QI+QAI+ICU+H+R+D+RD
  cumulative_infections <- rowSums(out[, index[ie_str %in% c(v_exp_str, v_asym_inf_str, v_inf_str,v_qe_str,v_qi_str,v_qai_str,
                                                                     "H", "ICU", "R","RD", "D")]])

  ## extract cumulative deaths by summing over deaths in all age groups and co-groups
  cumulative_deaths <- rowSums(out[, index[ie_str == "D"]])

  ## extract prevalent hospitalizations by adding people in ICU or H
  prevalent_hospitalizations <- rowSums(out[, index[ie_str %in% c("H", "ICU")]])
  # hospitalization for <30
  prevalent_hospitalizations1 <- rowSums(out[, index[ie_str %in% c("H") & ia %in% c(0:2)]])
  # hospitalization for 30-59
  prevalent_hospitalizations2 <- rowSums(out[, index[ie_str %in% c("H") & ia %in% c(3:5)]])
  # hospitalization for 60+
  prevalent_hospitalizations3 <- rowSums(out[, index[ie_str %in% c("H") & ia %in% c(6:8)]])
  
  ## proportion of 60 plus hospitalizaed
  tmp_60p <- rowSums(out[, index[ie_str %in% c("H", "ICU") & ia %in% c(6:8)]])
  prop_hosp_60p <- tmp_60p / prevalent_hospitalizations
  
  ## extract demand for ICU beds by adding people in the ICU in all age and co groups
  ICU_bed_demand <- rowSums(out[, index[ie_str == "ICU"]])

  ## calculate daily deaths given cumulative deaths
  daily_deaths <- c(0, diff(cumulative_deaths))
  
  ## extract cumulative home deaths by summing over deaths in all age groups and co-groups
  cumulative_home_deaths <- rowSums(out[, index[ie_str == "HD"]])
  ## Proportion of home deaths among 60+ who passed away
  prop_home_deaths_60p <- cumulative_home_deaths / rowSums(out[, index[ie_str == "D" & ia %in% c(6:8)]])
  ## Daily home death
  tmp_daily_death <- c(0, diff(cumulative_home_deaths))
  prop_hosp_60p_and_home_death <- (tmp_daily_death + tmp_60p) / (tmp_daily_death + prevalent_hospitalizations)
  
  ## extract cumulative hospitalizations by summing over deaths in all age groups and co-groups
  cumulative_hospitalizations <- rowSums(out[, index[ie_str == "CH"]])
  
  ## extract prevalent infections in I category alone
  prevalent_I_state <- rowSums(out[, index[ie_str %in% v_inf_str]])
  
  m_out <- cbind(prevalent_infections = prevalent_infections,
                 cumulative_infections = cumulative_infections, 
                 daily_deaths = daily_deaths,
                 cumulative_deaths = cumulative_deaths, 
                 ICU_bed_demand = ICU_bed_demand,
                 prevalent_hospitalizations = prevalent_hospitalizations,
                 prevalent_hospitalizations1 = prevalent_hospitalizations1, 
                 prevalent_hospitalizations2 = prevalent_hospitalizations2, 
                 prevalent_hospitalizations3 = prevalent_hospitalizations3, 
                 cumulative_home_deaths = cumulative_home_deaths, 
                 cumulative_hospitalizations = cumulative_hospitalizations,
                 prevalent_I_state = prevalent_I_state, 
                 prop_hosp_60p = prop_hosp_60p, 
                 prop_home_deaths_60p = prop_home_deaths_60p, 
                 prop_hosp_60p_and_home_death = prop_hosp_60p_and_home_death, 
                 out)
  
  out_data <- list()
  out_data$out <- m_out
  
  ## Arrange age group data
  if (chk) {
    name_data <- names(inc_data)
    ls_inc_data <- lapply(c(1:length(inc_data)), function(x) {
      
      tmp_name <- names(inc_data)[x]
      tmp_df <- inc_data[[tmp_name]]
      tmp_df <- data.table(tmp_df[-1, ])
      
      tmp_df[, `:=` (Time = (ceiling(Time) - 1))]
      tmp_df <- tmp_df[, lapply(.SD, sum), by = Time]
      tmp_df <- melt(tmp_df, 
                     id.vars = c("Time"), 
                     variable.name = "category", 
                     value.name = "incidence")
      tmp_df[, agp := tstrsplit(category, " ")[[3]]]
      tmp_df <- tmp_df[, list(inc = sum(incidence)), 
                       by = list(Time, agp)]
      
      return(tmp_df)  
    })
    names(ls_inc_data) <- name_data
    
    # incidence of I
    new_I <- ls_inc_data$new_I
    new_QI <- ls_inc_data$new_QI
    
    start_date <- set_start_date()
    
    incidence_data <- merge(new_I, new_QI, by = c("Time", "agp"), all = T)
    incidence_data[, `:=` (I = inc.x + inc.y, 
                           date = Time + start_date, 
                           age_gp = agp, 
                           inc.x = NULL, 
                           inc.y = NULL, 
                           agp = NULL, 
                           Time = NULL)]
    
    # incidence of AI
    new_AI <- ls_inc_data$new_AI
    new_QAI <- ls_inc_data$new_QAI
    
    df_AI <- merge(new_AI, new_QAI, by = c("Time", "agp"), all = T)
    df_AI[, `:=` (AI = inc.x + inc.y, 
                  date = Time + start_date, 
                  age_gp = agp, 
                  inc.x = NULL, 
                  inc.y = NULL, 
                  agp = NULL, 
                  Time = NULL)]
    incidence_data <- merge(incidence_data, df_AI, by = c("date", "age_gp"))
    
    new_H <- ls_inc_data$new_H 
    new_H[, `:=` (H = inc, 
                  date = Time + start_date, 
                  age_gp = agp, 
                  inc = NULL, 
                  agp = NULL, 
                  Time = NULL)]
    incidence_data <- merge(incidence_data, new_H, by = c("date", "age_gp"))
    
    new_ICU <- ls_inc_data$new_ICU
    new_ICU[, `:=` (ICU = inc, 
                    date = Time + start_date, 
                    age_gp = agp, 
                    inc = NULL, 
                    agp = NULL, 
                    Time = NULL)]
    incidence_data <- merge(incidence_data, new_ICU, by = c("date", "age_gp"))
    incidence_data[, `:=`(totI = I + AI, 
                          totH = H + ICU)]
    
    out_data$incidence_data <- incidence_data
    
    return(out_data)
  } else {
    return(out_data$out)
  }
}

#' @export
calculate_Rt <- function(dt, parms, end_date = NULL) {
  start_date <- set_start_date()
  if (is.null(end_date)) {
    end_t <- as.numeric(Sys.Date() - start_date) 
  } else {
    end_t <- as.numeric(end_date - start_date) 
  }
  # Rt estimation
  lm_Rt <- lm(log(dt[(end_t - 30) : end_t][["cumulative_infections"]]) ~ dt[(end_t - 30) : end_t][["Time"]])
  avg_exp_dur <- parms$n_exposed_states/(parms$exposed_transition_rate/parms$timestep)
  avg_inf_dur <- parms$n_infected_states/(parms$infected_transition_rate/parms$timestep)
  Rt_est <- (1 + lm_Rt$coefficients[[2]] * avg_inf_dur) * (1 + lm_Rt$coefficients[[2]] * avg_exp_dur)
  
  lm_R0 <- lm(log(dt[1 : 20][["cumulative_infections"]]) ~ dt[1 : 20][["Time"]])
  R0_est <- (1 + lm_R0$coefficients[[2]] * avg_inf_dur) * (1 + lm_R0$coefficients[[2]] * avg_exp_dur)
  return(c("R0" = R0_est, "Rt" = Rt_est))
}




#' @import data.table
#' @export
calculate_mdh_outcomes <- function(dt, parms, end_day_calibration) {
  R_est <- calculate_Rt(dt, parms)
  n_deaths <- round(dt[.N, "cumulative_deaths"][[1]], 0)
  n_deathsJulyEnd <- round(dt[date == as.Date("2020-07-31"), "cumulative_deaths"][[1]], 0)
  n_deathsAugEnd <- round(dt[date == as.Date("2020-08-31"), "cumulative_deaths"][[1]], 0)
  day_peak_deaths <- which.max(dt[["daily_deaths"]])
  date_peak_deaths <- dt[day_peak_deaths][["date"]]
  pct_deaths <- round(100 * n_deaths / parms$N, 2)
  day_peak_infections <- which.max(dt[["prevalent_infections"]])
  date_peak_infections <- dt[day_peak_infections][["date"]]
  peak_number_infections <- round(max(dt[["prevalent_infections"]]), 0)
  pct_infections <- round(100 * dt[.N, "cumulative_infections"][[1]] / parms$N, 2)
  day_peak_hospitalizations <- which.max(dt[["prevalent_hospitalizations"]])
  date_peak_hospitalizations <- dt[day_peak_hospitalizations][["date"]]
  peak_number_hospitalizations <- round(max(dt[["prevalent_hospitalizations"]]), 0)
  day_icu_cap_reached <- which(dt[["ICU_bed_demand"]] >= parms$n_icu_beds)[1]
  max_icu_demand <- round(max(dt[["ICU_bed_demand"]]), 0)
  prop_home_deaths_60p <- dt[end_day_calibration][["prop_home_deaths_60p"]]
  
  return(list("Rt_est" = R_est[["Rt"]], 
              "R0_est" = R_est[["R0"]], 
              "n_deaths" = n_deaths, "n_deathsJulyEnd" = n_deathsJulyEnd,
              "n_deathsAugEnd" = n_deathsAugEnd, 
              "day_peak_deaths" = day_peak_deaths, "date_peak_deaths" = date_peak_deaths, 
              "pct_deaths" = pct_deaths, "day_peak_infections" = day_peak_infections, 
              "date_peak_infections" = date_peak_infections, 
              "peak_number_infections" = peak_number_infections, 
              "pct_infections" = pct_infections, 
              "day_peak_hospitalizations" = day_peak_hospitalizations, 
              "date_peak_hospitalizations" = date_peak_hospitalizations, 
              "peak_number_hospitalizations" = peak_number_hospitalizations, 
              "day_icu_cap_reached" = day_icu_cap_reached, 
              "max_icu_demand" = max_icu_demand, 
              "prop_home_deaths_60p" = prop_home_deaths_60p))
}
