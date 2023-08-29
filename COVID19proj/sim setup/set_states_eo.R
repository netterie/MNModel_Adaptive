#' This is where we set up changes in policy, epi states, and hash_table of state indices
#' different simulations. 
#' The only input value required is the end_date. 
#' @param end_date The end day of the simulation. The end_date should be provided as a string, e.g., "2020-07-31"
#' @return The function return a list of predetermined parameters for simulations, including status_quo_strategy, 
#' epi_states_index, hash_table
set_states_eo <- function(end_date, end_date_sip = as.Date("2020-05-18")) {
  start_date <- COVID19pack::set_start_date()
  
  status_quo_strategy <- list(start_time_social_distancing = as.numeric(as.Date("2020-03-23") - start_date), 
                              end_time_social_distancing = as.numeric(as.Date("2020-03-28") - start_date),  
                              start_time_sip = as.numeric(as.Date("2020-03-28") - start_date), 
                              end_time_sip = as.numeric(end_date_sip - start_date), 
                              start_time_60plus_distancing = as.numeric(end_date_sip - start_date), 
                              end_time_60plus_distancing = as.numeric(as.Date(end_date) - start_date), 
                              start_time_behavior_change = as.numeric(end_date_sip - start_date),
                              end_time_behavior_change = as.numeric(as.Date(end_date) - start_date)) 
  
  #### Set up the epi states and indices
  epi_states_index <- gen_states_index(age_groups = c(0, 10, 20, 30, 40, 50, 60, 70, 80),
                                       comorbidity_groups = c("c0", "c1"), 
                                       n_infected_states = 3, 
                                       n_exposed_states = 2)
  
  ## create hash_table
  hash_table <- create_hash_table(epi_states_index)
  
  ls_return <- list(start_date = start_date, status_quo_strategy = status_quo_strategy, 
                    epi_states_index = epi_states_index, hash_table = hash_table)
  
  
  return(ls_return)
}
