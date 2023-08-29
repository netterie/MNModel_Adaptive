# This is a script to replicate model calibration


library(EasyABC)
library(COVID19pack)
library(readxl)

source("COVID19proj/sim setup/set_states_eo.R")

### ABC process
sim_output <- function(x, 
                       calibrate = TRUE) {
  
  
  start_date <- set_start_date()
  dat_date <- as.Date("2020-08-09")
  end_date <- as.Date("2021-03-22")
  cal_date <- as.Date("2020-12-01")
  end_day <- as.numeric(end_date - start_date)
  dat_day <- as.numeric(dat_date - start_date)
  cal_day <- as.numeric(cal_date - start_date)
  tmp_ls <- set_states_eo(end_date, as.Date("2020-05-18"))
  list2env(tmp_ls, envir = .GlobalEnv)
  
  parms <- COVID19pack:::gen_abc_param23(x[1:17], 
                                         epi_states_index = epi_states_index, 
                                         hash_table = hash_table)
  parms$trigger_strategy <- list(trigger_thres = c(icu = 2200 * 0.1, case2x = 7, 
                                                   hosp = 100000/100000, prev = 0.05), 
                                 trigger = "hosp", 
                                 days_past_peak = 14, #days
                                 grace_period = 7, #days
                                 min_on_day = 14, #days
                                 on_contact_red = 0.8,
                                 off_contact_red = x[18],
                                 on_contact_red_60p = 0.8,
                                 off_contact_red_60p = x[18])
  
  names_sq_st <- names(status_quo_strategy)
  for (i in names_sq_st) {
    parms[[i]] <- status_quo_strategy[[i]]
  }
  
  parms$str_peak_type <- "hospitalizations"
  parms$end_time_behavior_change <- as.numeric(as.Date("2021-03-22") - start_date)
  parms$end_time_60plus_distancing <- as.numeric(as.Date("2020-09-01") - start_date)
  
  times <- seq(1, 365, by = parms$timestep)
  
  m_out_raw <- solve_model(times = times,
                           func = covid_19_model_function,
                           parms = parms, 
                           return_full = 0, 
                           return_incidence = TRUE)
  res <- process_output(m_out_raw, parms)
  sim_data <- data.table(res$out)
  sim_data[, date := (Time + start_date)]
  tix <- seq(7, dat_day, 7)
  tix2 <- seq(dat_day, cal_day, 14)
  
  sim_out <- round(c(sim_data$cumulative_deaths[tix], 
                     sim_data$prevalent_hospitalizations1[tix], 
                     sim_data$prevalent_hospitalizations2[tix], 
                     sim_data$prevalent_hospitalizations3[tix],
                     sim_data$prevalent_hospitalizations[tix2]))
  return(sim_out)
}

# Read in Calibration Targets
target_ls <- readRDS("COVID19proj/data/Calib_Targets.rds")

parm_bd_df <- COVID19pack:::spec_parm23()
parm_bd_df[nrow(parm_bd_df)+1,] <- list("off_cr", lb = 0.01, ub = 0.99)

## Set end day
dat_day <- as.numeric(as.Date("2020-08-09") - start_date)
cal_day <- as.numeric(as.Date("2020-12-01") - start_date)

## generate targets
tix <- seq(7, dat_day, 7)
tix2 <- seq(dat_day, cal_day, 14)
targets <- c(target_ls$actual_deaths[tix], 
             target_ls$actual_hospitalizations1[tix], 
             target_ls$actual_hospitalizations2[tix], 
             target_ls$actual_hospitalizations3[tix],
             rep(377, length(tix2)))



priors <- lapply(c(1:nrow(parm_bd_df)), function(x) {
  c("unif", parm_bd_df[x, "lb"], parm_bd_df[x, "ub"])
})


nsim <- 1000

aa <- Sys.time()
abs_res <- ABC_sequential(method = "Lenormand",
                          model = sim_output,
                          prior = priors,
                          nb_simul = nsim,
                          summary_stat_target = targets,
                          p_acc_min = 0.05, 
                          # dist_weights = target_wt, 
                          progress_bar = F) # , 
# n_cluster = ncore,
# use_seed = TRUE)
print(Sys.time() - aa)

# Save the calibrated parameters 
parm_df <- data.frame(abs_res$param)
parm_df$weights <- abs_res$weights

# Save posterior parameters
saveRDS(parm_df, file = "COVID19proj/data/abc_posterior_samples.rds")


