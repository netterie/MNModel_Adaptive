# ****************************************************************************
# Reduction of Scenarios.R to only two scenarios and "nsamp" posterior samples,
# using tidied file structure and path references
# ****************************************************************************

# Assumes you are in COVID19_Trigger.Rproj
# set runmodel=FALSE to do the setup, skip running the model, and return "parms"
# uncomment the line at the end of this file to run scenarios A and B

test_scenarios = function(nsamp=1, runmodel=FALSE) {
  
  # ----------------------------------------------------------------------------
  # SETUP
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
  
  # Source helper functions in the setup folder
  setup_files = as.character(file.path(dir_setup, list.files(dir_setup)))
  sapply(setup_files, source)
  
  # Load posterior samples 
  x<-readRDS(file.path(dir_data, "abc_posterior_samples.rds"))
  
  # Prepare dates
  # ---Start
  start_date <- set_start_date()
  # ---Dat
  dat_date <- as.Date("2020-08-09")
  dat_day <- as.numeric(dat_date - start_date)
  # ---End
  end_date <- as.Date("2021-03-22")
  end_day <- as.numeric(end_date - start_date)
  
  # ----------------------------------------------------------------------------
  # PREPARE ENVIRONMENT
  # ----------------------------------------------------------------------------
  
  # Load ..defaults? into environment
  tmp_ls <- set_states_eo(end_date)
  list2env(tmp_ls, envir = .GlobalEnv)
  
  # Prepare results
  results_A <- data.frame(A1=c(1:365))
  results_B <- data.frame(A1=c(1:365))
  
  NPI_A <- data.frame(A1=c(1:365))
  NPI_B <- data.frame(A1=c(1:365))
  
  # ----------------------------------------------------------------------------
  # LOOP THROUGH POSTERIOR SAMPLES AND RUN MODEL
  # ----------------------------------------------------------------------------
  
  # x contains the posterior samples for the estimated parameters
  # Each row of x is one sampled parameter set
  # Each scenario is run using each parameter set
  # x is generated by COVID19proj/code/CalibABC.R
  for(j in c(1:nsamp)){
    
    # Use the package to assign parameter set j to the correct 
    # parameter values in parms, the default settings.
    # Requires a function because the 
    # parameters in x need to be accessed as-is or multiplied in 
    # various combinations and structures (vectors vs matrices)
    parms <- COVID19pack:::gen_abc_param23(x[j,1:17], 
                                           epi_states_index = epi_states_index, 
                                           hash_table = hash_table)
    names_sq_st <- names(status_quo_strategy)
    # Establish status quo strategy in parms
    for (i in names_sq_st) {
      parms[[i]] <- status_quo_strategy[[i]]
    }
    parms$str_peak_type <- "hospitalizations"
    parms$end_time_behavior_change <- as.numeric(as.Date("2021-03-22") - start_date)
    parms$end_time_60plus_distancing <- as.numeric(as.Date("2020-09-01") - start_date)
    
    parms_A<- parms_B <- parms
    
    times <- seq(1, 365, by = parms$timestep)
    
    ## Scenario
    # Steady state (Scenario A):
    parms_A$trigger_strategy <- list(trigger_thres = c(icu = 2200 * 0.1, case2x = 7, 
                                                       hosp = 10000/100000, prev = 0.05), 
                                     trigger = "hosp", 
                                     days_past_peak = 14, #days
                                     grace_period = 7, #days
                                     min_on_day = 14, #days
                                     on_contact_red = 0.8,
                                     off_contact_red = x[j,18],
                                     on_contact_red_60p = 0.8,
                                     off_contact_red_60p = x[j,18])
    
    # COVID fatigue (Scenario B):
    parms_B$trigger_strategy <- list(trigger_thres = c(icu = 2200 * 0.1, case2x = 7, 
                                                       hosp = 10000/100000, prev = 0.05), 
                                     trigger = "hosp", 
                                     days_past_peak = 14, #days
                                     grace_period = 7, #days
                                     min_on_day = 14, #days
                                     on_contact_red = 0.8,
                                     off_contact_red = x[j,18] - 0.1,
                                     on_contact_red_60p = 0.8,
                                     off_contact_red_60p = x[j,18] - 0.1)
    
    if (runmodel) {
      
      # Scenario A
      # Solve using custom built solver
      m_out_raw_A <- solve_model(times = times,
                                 func = covid_19_model_function,
                                 parms = parms_A, 
                                 return_full = 0, 
                                 return_incidence = TRUE)
      
      res_A <- process_output(m_out_raw_A, parms_A)
      sim_data_A <- data.table(res_A$out)
      sim_data_A[, date := (Time + start_date)]
      results_A[,j+1] <- sim_data_A$prevalent_hospitalizations
      #### calculate mdh outcomes
      sim_data_A[, NPI := ifelse(sd==1,"sd",NA)]
      sim_data_A[, NPI := ifelse(sip==1,"SAH",NPI)]
      sim_data_A[, NPI := ifelse(bec==1,"Post-SAH",NPI)]
      sim_data_A[date>"2020-09-01"&NPI=="SAH", NPI := "trigger_on"]
      sim_data_A[date>"2020-09-01"&NPI=="Post-SAH", NPI := "trigger_off"]
      sim_data_A$NPI <- factor(sim_data_A$NPI, 
                               levels=c("sd", "SAH", "Post-SAH", "trigger_off", "trigger_on", "NA"), 
                               labels=c("Period 1", "Period 2", "Period 3", "50% effective contact reduction", 
                                        "80% effective contact reduction", "NA"))
      NPI_A[,j+1] <- sim_data_A$NPI 
      
      # Scenario B
      # Solve using custom built solver
      m_out_raw_B <- solve_model(times = times,
                                 func = covid_19_model_function,
                                 parms = parms_B, 
                                 return_full = 0, 
                                 return_incidence = TRUE)
      res_B <- process_output(m_out_raw_B, parms_B)
      sim_data_B <- data.table(res_B$out)
      sim_data_B[, date := (Time + start_date)]
      results_B[,j+1] <- sim_data_B$prevalent_hospitalizations
      #### calculate mdh outcomes
      sim_data_B[, NPI := ifelse(sd==1,"sd",NA)]
      sim_data_B[, NPI := ifelse(sip==1,"SAH",NPI)]
      sim_data_B[, NPI := ifelse(bec==1,"Post-SAH",NPI)]
      sim_data_B[date>"2020-09-01"&NPI=="SAH", NPI := "trigger_on"]
      sim_data_B[date>"2020-09-01"&NPI=="Post-SAH", NPI := "trigger_off"]
      sim_data_B$NPI <- factor(sim_data_B$NPI, 
                               levels=c("sd", "SAH", "Post-SAH", "trigger_off", "trigger_on", "NA"), 
                               labels=c("Period 1", "Period 2", "Period 3", "50% effective contact reduction", 
                                        "80% effective contact reduction", "NA"))
      NPI_B[,j+1] <- sim_data_B$NPI 
    } #end if (runmodel)
    
    print(j)
  }
  
  # ----------------------------------------------------------------------------
  # SAVE RESULTS
  # ----------------------------------------------------------------------------
  if (runmodel) {
    for(k in c("A", "B")){
      
      saveRDS(get(paste0("results_", k))[,-1], file = paste0("data/CalRes_", k ,".rds"))
      saveRDS(get(paste0("NPI_", k))[,-1], file = paste0("data/Trigger_", k ,".rds"))
      
    } #end for(k...)
    
    return()
    
  } else return(parms)

} #end test_scenarios


# test_scenarios(nsamp=1, runmodel=TRUE)
