# This file is to generate base case results dataset for scenario A-F

rm(list = ls())

# Loading COVID-19 package
library(COVID19pack) # make sure you have installed the COVID19pack package
library(ggplot2)
library(reshape2)

#### Set working directory (Please ensure your current working directory is under MNCOVID19)
# setwd(paste0(getwd(), "/COVID19proj"))

x<-readRDS("data/abc_posterior_samples.rds")

source("sim setup/set_states_eo.R") # helper function for running the model

start_date <- set_start_date()
dat_date <- as.Date("2020-08-09")
end_date <- as.Date("2021-03-22")
end_day <- as.numeric(end_date - start_date)
dat_day <- as.numeric(dat_date - start_date)
tmp_ls <- set_states_eo(end_date)
list2env(tmp_ls, envir = .GlobalEnv)

results_A <- data.frame(A1=c(1:365))
results_B <- data.frame(A1=c(1:365))
results_C <- data.frame(A1=c(1:365))
results_D <- data.frame(A1=c(1:365))
results_E <- data.frame(A1=c(1:365))
results_F <- data.frame(A1=c(1:365))

NPI_A <- data.frame(A1=c(1:365))
NPI_B <- data.frame(A1=c(1:365))
NPI_C <- data.frame(A1=c(1:365))
NPI_D <- data.frame(A1=c(1:365))
NPI_E <- data.frame(A1=c(1:365))
NPI_F <- data.frame(A1=c(1:365))



for(j in c(1:nrow(x))){
  
  parms <- COVID19pack:::gen_abc_param23(x[j,1:17], 
                                         epi_states_index = epi_states_index, 
                                         hash_table = hash_table)
  names_sq_st <- names(status_quo_strategy)
  for (i in names_sq_st) {
    parms[[i]] <- status_quo_strategy[[i]]
  }
  parms$str_peak_type <- "hospitalizations"
  parms$end_time_behavior_change <- as.numeric(as.Date("2021-03-22") - start_date)
  parms$end_time_60plus_distancing <- as.numeric(as.Date("2020-09-01") - start_date)
  
  parms_A<- parms_B <-parms_C <- parms_D <- parms_E<- parms_F<- parms
  
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
  
  # Trigger reference (Scenario C):
  parms_C$trigger_strategy <- list(trigger_thres = c(icu = 2200 * 0.1, case2x = 7, 
                                                     hosp = 8/100000, prev = 0.05), 
                                   trigger = "hosp", 
                                   days_past_peak = 14, #days
                                   grace_period = 7, #days
                                   min_on_day = 14, #days
                                   on_contact_red = 0.8,
                                   off_contact_red = x[j,18] - 0.1,
                                   on_contact_red_60p = 0.8,
                                   off_contact_red_60p = x[j,18] - 0.1)
  
  # Slower adjustment (Scenario D):
  parms_D$trigger_strategy <- list(trigger_thres = c(icu = 2200 * 0.1, case2x = 7, 
                                                     hosp = 8/100000, prev = 0.05), 
                                   trigger = "hosp", 
                                   days_past_peak = 14, #days
                                   grace_period = 14, #days
                                   min_on_day = 14, #days
                                   on_contact_red = 0.8,
                                   off_contact_red = x[j,18] - 0.1,
                                   on_contact_red_60p = 0.8,
                                   off_contact_red_60p = x[j,18] - 0.1)
  
  # Stronger response (Scenario E):
  parms_E$trigger_strategy <- list(trigger_thres = c(icu = 2200 * 0.1, case2x = 7, 
                                                     hosp = 8/100000, prev = 0.05), 
                                   trigger = "hosp", 
                                   days_past_peak = 14, #days
                                   grace_period = 7, #days
                                   min_on_day = 14, #days
                                   on_contact_red = 0.9,
                                   off_contact_red = x[j,18] - 0.1,
                                   on_contact_red_60p = 0.9,
                                   off_contact_red_60p = x[j,18] - 0.1)
  
  # Higher threshold (Scenario F):
  parms_F$trigger_strategy <- list(trigger_thres = c(icu = 2200 * 0.1, case2x = 7, 
                                                     hosp = 10/100000, prev = 0.05), 
                                   trigger = "hosp", 
                                   days_past_peak = 14, #days
                                   grace_period = 7, #days
                                   min_on_day = 14, #days
                                   on_contact_red = 0.8,
                                   off_contact_red = x[j,18] - 0.1,
                                   on_contact_red_60p = 0.8,
                                   off_contact_red_60p = x[j,18] - 0.1)
  
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
  
  # Scenario C
  # Solve using custom built solver
  m_out_raw_C <- solve_model(times = times,
                             func = covid_19_model_function,
                             parms = parms_C, 
                             return_full = 0, 
                             return_incidence = TRUE)
  res_C <- process_output(m_out_raw_C, parms_C)
  sim_data_C <- data.table(res_C$out)
  sim_data_C[, date := (Time + start_date)]
  results_C[,j+1] <- sim_data_C$prevalent_hospitalizations
  #### calculate mdh outcomes
  sim_data_C[, NPI := ifelse(sd==1,"sd",NA)]
  sim_data_C[, NPI := ifelse(sip==1,"SAH",NPI)]
  sim_data_C[, NPI := ifelse(bec==1,"Post-SAH",NPI)]
  sim_data_C[date>"2020-09-01"&NPI=="SAH", NPI := "trigger_on"]
  sim_data_C[date>"2020-09-01"&NPI=="Post-SAH", NPI := "trigger_off"]
  sim_data_C$NPI <- factor(sim_data_C$NPI, 
                           levels=c("sd", "SAH", "Post-SAH", "trigger_off", "trigger_on", "NA"), 
                           labels=c("Period 1", "Period 2", "Period 3", "50% effective contact reduction", 
                                    "80% effective contact reduction", "NA"))
  NPI_C[,j+1] <- sim_data_C$NPI 
  
  # Scenario D
  # Solve using custom built solver
  m_out_raw_D <- solve_model(times = times,
                             func = covid_19_model_function,
                             parms = parms_D, 
                             return_full = 0, 
                             return_incidence = TRUE)
  res_D <- process_output(m_out_raw_D, parms_D)
  sim_data_D <- data.table(res_D$out)
  sim_data_D[, date := (Time + start_date)]
  results_D[,j+1] <- sim_data_D$prevalent_hospitalizations
  #### calculate mdh outcomes
  sim_data_D[, NPI := ifelse(sd==1,"sd",NA)]
  sim_data_D[, NPI := ifelse(sip==1,"SAH",NPI)]
  sim_data_D[, NPI := ifelse(bec==1,"Post-SAH",NPI)]
  sim_data_D[date>"2020-09-01"&NPI=="SAH", NPI := "trigger_on"]
  sim_data_D[date>"2020-09-01"&NPI=="Post-SAH", NPI := "trigger_off"]
  sim_data_D$NPI <- factor(sim_data_D$NPI, 
                           levels=c("sd", "SAH", "Post-SAH", "trigger_off", "trigger_on", "NA"), 
                           labels=c("Period 1", "Period 2", "Period 3", "50% effective contact reduction", 
                                    "80% effective contact reduction", "NA"))
  NPI_D[,j+1] <- sim_data_D$NPI 
  
  # Scenario E
  # Solve using custom built solver
  m_out_raw_E <- solve_model(times = times,
                             func = covid_19_model_function,
                             parms = parms_E, 
                             return_full = 0, 
                             return_incidence = TRUE)
  res_E <- process_output(m_out_raw_E, parms_E)
  sim_data_E <- data.table(res_E$out)
  sim_data_E[, date := (Time + start_date)]
  results_E[,j+1] <- sim_data_E$prevalent_hospitalizations
  #### calculate mdh outcomes
  sim_data_E[, NPI := ifelse(sd==1,"sd",NA)]
  sim_data_E[, NPI := ifelse(sip==1,"SAH",NPI)]
  sim_data_E[, NPI := ifelse(bec==1,"Post-SAH",NPI)]
  sim_data_E[date>"2020-09-01"&NPI=="SAH", NPI := "trigger_on"]
  sim_data_E[date>"2020-09-01"&NPI=="Post-SAH", NPI := "trigger_off"]
  sim_data_E$NPI <- factor(sim_data_E$NPI, 
                           levels=c("sd", "SAH", "Post-SAH", "trigger_off", "trigger_on", "NA"), 
                           labels=c("Period 1", "Period 2", "Period 3", "50% effective contact reduction", 
                                    "80% effective contact reduction", "NA"))
  NPI_E[,j+1] <- sim_data_E$NPI 
  
  # Scenario F
  # Solve using custom built solver
  m_out_raw_F <- solve_model(times = times,
                             func = covid_19_model_function,
                             parms = parms_F, 
                             return_full = 0, 
                             return_incidence = TRUE)
  res_F <- process_output(m_out_raw_F, parms_F)
  sim_data_F <- data.table(res_F$out)
  sim_data_F[, date := (Time + start_date)]
  results_F[,j+1] <- sim_data_F$prevalent_hospitalizations
  #### calculate mdh outcomes
  sim_data_F[, NPI := ifelse(sd==1,"sd",NA)]
  sim_data_F[, NPI := ifelse(sip==1,"SAH",NPI)]
  sim_data_F[, NPI := ifelse(bec==1,"Post-SAH",NPI)]
  sim_data_F[date>"2020-09-01"&NPI=="SAH", NPI := "trigger_on"]
  sim_data_F[date>"2020-09-01"&NPI=="Post-SAH", NPI := "trigger_off"]
  sim_data_F$NPI <- factor(sim_data_F$NPI, 
                           levels=c("sd", "SAH", "Post-SAH", "trigger_off", "trigger_on", "NA"), 
                           labels=c("Period 1", "Period 2", "Period 3", "50% effective contact reduction", 
                                    "80% effective contact reduction", "NA"))
  NPI_F[,j+1] <- sim_data_F$NPI 
  
  
  print(j)
}

for(k in c("A", "B", "C", "D", "E", "F")){
  
  saveRDS(get(paste0("results_", k))[,-1], file = paste0("data/CalRes_", k ,".rds"))
  saveRDS(get(paste0("NPI_", k))[,-1], file = paste0("data/Trigger_", k ,".rds"))
  
}
