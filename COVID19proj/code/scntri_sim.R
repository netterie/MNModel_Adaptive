# This file includes scntri_sim function, to calculate results for sensitivity analysis of 
# triggering threshold, adjustment time and NPI effectiveness (shown in Figure S3)

scntri_sim <- function(name, 
                       eff = 0.8,
                       fatigue = 0.1,
                       thresh = 8, 
                       grace_period = 7, 
                       beta_change = FALSE,
                       delta_beta_after = 1){
  
  x<-readRDS("data/abc_samples_RR22.rds")
  parm_bd_df <- COVID19pack:::spec_parm23()
  
  source("sim setup/set_states_eo.R") # helper function for running the model
  
  start_date <- set_start_date()
  dat_date <- as.Date("2020-08-09")
  end_date <- as.Date("2021-03-22")
  end_day <- as.numeric(end_date - start_date)
  dat_day <- as.numeric(dat_date - start_date)
  tmp_ls <- set_states_eo(end_date)
  list2env(tmp_ls, envir = .GlobalEnv)
  
  results <- data.frame(A1=c(1:365))

  
  NPI <- data.frame(A1=c(1:365))


  
  for(j in c(1:nrow(x))){
    
    parms <- COVID19pack:::gen_abc_param23(x[j,1:17], 
                                           epi_states_index = epi_states_index, 
                                           hash_table = hash_table,
                                           beta_change = beta_change,
                                           beta_after = x[j,6] * delta_beta_after)
    names_sq_st <- names(status_quo_strategy)
    for (i in names_sq_st) {
      parms[[i]] <- status_quo_strategy[[i]]
    }
    parms$str_peak_type <- "hospitalizations"
    parms$end_time_behavior_change <- as.numeric(as.Date("2021-03-22") - start_date)
    parms$end_time_60plus_distancing <- as.numeric(as.Date("2020-09-01") - start_date)
    
    times <- seq(1, 365, by = parms$timestep)
    
    ## Scenario
    # Trigger reference (Scenario C):
    parms$trigger_strategy <- list(trigger_thres = c(icu = 2200 * 0.1, case2x = 7, 
                                                     hosp = thresh/100000, prev = 0.05), 
                                     trigger = "hosp", 
                                     days_past_peak = 14, #days
                                     grace_period = grace_period, #days
                                     min_on_day = 14, #days
                                     on_contact_red = eff,
                                     off_contact_red = x[j,18] - fatigue,
                                     on_contact_red_60p = eff,
                                     off_contact_red_60p = x[j,18] - fatigue)
    
    # Scenario C
    # Solve using custom built solver
    m_out_raw  <-  solve_model(times = times,
                               func = covid_19_model_function,
                               parms = parms, 
                               return_full = 0, 
                               return_incidence = TRUE)
    res <- process_output(m_out_raw, parms)
    sim_data <- data.table(res$out)
    sim_data[, date := (Time + start_date)]
    results[,j+1] <- sim_data$prevalent_hospitalizations
    #### calculate mdh outcomes
    sim_data[, NPI := ifelse(sd==1,"sd",NA)]
    sim_data[, NPI := ifelse(sip==1,"SAH",NPI)]
    sim_data[, NPI := ifelse(bec==1,"Post-SAH",NPI)]
    sim_data[date>"2020-09-01"&NPI=="SAH", NPI := "trigger_on"]
    sim_data[date>"2020-09-01"&NPI=="Post-SAH", NPI := "trigger_off"]
    sim_data$NPI <- factor(sim_data$NPI, 
                             levels=c("sd", "SAH", "Post-SAH", "trigger_off", "trigger_on", "NA"), 
                             labels=c("Period 1", "Period 2", "Period 3", "50% effective contact reduction", 
                                      "80% effective contact reduction", "NA"))
    NPI[,j+1] <- sim_data$NPI 
    print(j)
  }

    saveRDS(get("results")[,-1], file = paste0("data/sa/CalRes_", name,".rds"))
    saveRDS(get("NPI")[,-1], file = paste0("data/sa/Trigger_", name,".rds"))

}