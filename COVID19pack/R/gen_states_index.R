#' Generate epi states and their indecies
#'
#' @param age_groups a vector of age group
#' @param comorbidity_groups a string vector of different comorbidity groups
#' @param n_exposed_states number of exposed states
#' @param n_infected_states number of infected states
#'
#' @return
#' a list of fixed population characteristics for the simulation
#' 
#' @export
gen_states_index <- function(age_groups = c(0, 10, 20, 30, 40, 50, 60, 70, 80),
                             comorbidity_groups = c("c0", "c1"), 
                             n_exposed_states, 
                             n_infected_states) {
  #total population#
  N <- sum(population["pop"])
  
  n_age_groups <- length(age_groups)
  
  N_by_age <- rep(0, n_age_groups)
  
  for (k in 1:(length(age_groups) - 1)) {
    N_by_age[k] <- sum(population[(age_groups[k]):(age_groups[k + 1] - 1), 2])
  }
  N_by_age[n_age_groups] <- sum(population[(age_groups[n_age_groups]):nrow(population), 2])
  age_prop <- N_by_age / N
  
  ncg <- length(comorbidity_groups)
  
  v_exp_str <- paste0("E", 1:n_exposed_states)
  v_qe_str <- paste0("QE",1:n_exposed_states)
  v_inf_str <- paste0("I", 1:n_infected_states)
  v_asym_inf_str <- paste0("AI", 1:n_infected_states)
  v_qi_str <- paste0("QI", 1:n_infected_states)
  v_qai_str <- paste0("QAI", 1:n_infected_states)
  epi_groups <- c("S", v_exp_str, v_asym_inf_str, v_inf_str,
                  "QS",v_qe_str, v_qai_str, v_qi_str, 
                  "H", "ICU", "R","RD","D", "HD", "CH")
  n_epi_groups <- length(epi_groups)
  
  # JKB NOTE: below and elsewhere, use of individual vectors and lists 
  # makes it difficult to quickly understand which vectors
  # have the same dimensions and correspond to the same conceptual value,
  # whereas a data frame makes that obvious
  
  ## Names for output objects
  # JKB NOTE: This approach puts the whole indexing process at the mercy of 
  # the "outer" command
  nam <- c(outer(outer(comorbidity_groups, epi_groups, FUN = "paste"), age_groups , FUN = "paste"))
  
  ## Create index list
  # JKB NOTE: I believe this process of getting numeric indices is more complicated
  # than it needed to be because it had to match what was done by the "outer" command
  # ls_ix a list of
  # 1. ie = vector of the numeric indices (1:25) of the epi group for each of the 850 compartments
  # 2. ic = vector of the numeric indices (0:1) of the comorbidity group for each of the 850 compartments
  # 3. ia = vector of the numeric indices (0:16) of the age group for each of the 850 compartments
  # 4. ie_str = vector of the alphanumeric ids of the epi group for each of the 850 compartments (25 unique vals)
  # 4. index = vector of the numeric indices (1:850) of each compartment, but NOT in numeric order because
  #    of the formula approach used (see get_ind function below, "index <- eg * ncg + cg + ag * neg * ncg - (ncg - 1)"),
  #    probably because the 
  ls_ix <- get_ind(ncg = ncg, neg = n_epi_groups, nag = n_age_groups, epi_groups)
  
  epi_groups_ls <- list(v_exp_str = v_exp_str,
                        v_qe_str = v_qe_str,
                        v_inf_str = v_inf_str,
                        v_asym_inf_str = v_asym_inf_str,
                        v_qi_str = v_qi_str,
                        v_qai_str = v_qai_str)
  
  return(list(N = N, 
              n_age_groups = n_age_groups, 
              age_groups = age_groups, 
              N_by_age = N_by_age, 
              age_prop = age_prop, 
              ncg = ncg, 
              n_exposed_states = n_exposed_states, 
              n_infected_states = n_infected_states, 
              n_epi_groups = n_epi_groups, 
              epi_groups_ls = epi_groups_ls, 
              nam = nam, 
              ls_ix = ls_ix))
}

#' @export
create_hash_table <- function(epi_states_index) {
  
  # Extract numbers of states, strata, infected and exposed states
  nag <- epi_states_index$n_age_groups
  ncg <- epi_states_index$ncg
  nes <- epi_states_index$n_exposed_states
  nis <- epi_states_index$n_infected_states
  
  # Numeric indexes for states, strata, and each compartment. Plus ie_str (string ID)
  ie <- epi_states_index$ls_ix$ie
  ic <- epi_states_index$ls_ix$ic
  ia <- epi_states_index$ls_ix$ia
  ie_str <- epi_states_index$ls_ix$ie_str
  index <- epi_states_index$ls_ix$index

  for (i in c(1:length(epi_states_index$epi_groups_ls))) {
    assign(names(epi_states_index$epi_groups_ls)[i], epi_states_index$epi_groups_ls[[i]])
  }
    
  e <- new.env()
  
  #### indices for difference vector "d" in covid_19_model_function_v3
  e$v_S_ind <- index[ie_str == "S"]
  e$v_QS_ind <- index[ie_str == "QS"]
  e$v_E_ind <- index[ie_str %in% v_exp_str]
  e$v_E_ind1 <- e$v_E_ind[seq(1, length(e$v_E_ind), nes)]
  e$v_E_ind_rm_last <- e$v_E_ind[setdiff(c(1:length(e$v_E_ind)), seq(nes, length(e$v_E_ind), nes))]
  e$v_E_ind_rm_first <- e$v_E_ind[setdiff(c(1:length(e$v_E_ind)), seq(1, length(e$v_E_ind), nes))]
  e$v_E_ind_last <- e$v_E_ind[seq(nes, length(e$v_E_ind), nes)]
  e$v_AI_ind <- index[ie_str %in% v_asym_inf_str]
  e$v_AI_ind1 <- e$v_AI_ind[seq(1, length(e$v_AI_ind), nis)]
  e$v_AI_ind_rm_last <- e$v_AI_ind[setdiff(c(1:length(e$v_AI_ind)), seq(nis, length(e$v_AI_ind), nis))]
  e$v_AI_ind_rm_first <- e$v_AI_ind[setdiff(c(1:length(e$v_AI_ind)), seq(1, length(e$v_AI_ind), nis))]
  e$v_AI_ind_last <- e$v_AI_ind[seq(nis, length(e$v_AI_ind), nis)]
  e$v_I_ind <- index[ie_str %in% v_inf_str]
  e$v_I_ind1 <- e$v_I_ind[seq(1, length(e$v_I_ind), nis)]
  e$v_I_ind_rm_last <- e$v_I_ind[setdiff(c(1:length(e$v_I_ind)), seq(nis, length(e$v_I_ind), nis))]
  e$v_I_ind_rm_first <- e$v_I_ind[setdiff(c(1:length(e$v_I_ind)), seq(1, length(e$v_I_ind), nis))]
  e$v_I_ind_last <- e$v_I_ind[seq(nis, length(e$v_I_ind), nis)]
  e$v_QS_ind <- index[ie_str == "QS"]
  e$v_QE_ind <- index[ie_str %in% v_qe_str]
  e$v_QE_ind1 <- e$v_QE_ind[seq(1, length(e$v_QE_ind), nes)]
  e$v_QE_ind_rm_last <- e$v_QE_ind[setdiff(c(1:length(e$v_QE_ind)), seq(nes, length(e$v_QE_ind), nes))]
  e$v_QE_ind_rm_first <- e$v_QE_ind[setdiff(c(1:length(e$v_QE_ind)), seq(1, length(e$v_QE_ind), nes))]
  e$v_QE_ind_last <- e$v_QE_ind[seq(nes, length(e$v_QE_ind), nes)]
  e$v_QAI_ind <- index[ie_str %in% v_qai_str]
  e$v_QAI_ind1 <- e$v_QAI_ind[seq(1, length(e$v_QAI_ind), nis)]
  e$v_QAI_ind_rm_last <- e$v_QAI_ind[setdiff(c(1:length(e$v_QAI_ind)), seq(nis, length(e$v_QAI_ind), nis))]
  e$v_QAI_ind_rm_first <- e$v_QAI_ind[setdiff(c(1:length(e$v_QAI_ind)), seq(1, length(e$v_QAI_ind), nis))]
  e$v_QAI_ind_last <- e$v_QAI_ind[seq(nis, length(e$v_QAI_ind), nis)]
  e$v_QI_ind <- index[ie_str %in% v_qi_str]
  e$v_QI_ind1 <- e$v_QI_ind[seq(1, length(e$v_QI_ind), nis)]
  e$v_QI_ind_rm_last <- e$v_QI_ind[setdiff(c(1:length(e$v_QI_ind)), seq(nis, length(e$v_QI_ind), nis))]
  e$v_QI_ind_rm_first <- e$v_QI_ind[setdiff(c(1:length(e$v_QI_ind)), seq(1, length(e$v_QI_ind), nis))]
  e$v_QI_ind_last <-  e$v_QI_ind[seq(nis, length(e$v_QI_ind), nis)]
  e$v_H_ind <- index[ie_str == "H"]
  e$v_ICU_ind <- index[ie_str == "ICU"]
  e$v_R_ind <- index[ie_str == "R"]
  e$v_RD_ind <- index[ie_str == "RD"]
  e$v_D_ind <- index[ie_str == "D"]
  e$v_HD_ind <- index[ie_str == "HD"]
  e$v_CH_ind <- index[ie_str == "CH"]
  
  #### indices for lambda calculation in lambda.R file
  e$v_ind_inf_by_age <- index[ie_str %in% c(v_asym_inf_str, v_inf_str)]
  e$v_ind_D_by_age <- index[ie_str == "D"]
  e$v_ind_quar_inf_by_age <- index[ie_str %in% c(v_qai_str, v_qi_str)]
  e$v_ind_quar_noninf_by_age <- index[ie_str %in% c("QS",v_qe_str,"QR")]
  
  #### indices for increase in testing (originally in covid_19_model_function())
  e$s_index <- index[ie_str == "S"]
  e$exp_index <- index[ie_str %in% v_exp_str]
  e$inf_index <- index[ie_str %in% v_inf_str]
  e$asym_index <- index[ie_str %in% v_asym_inf_str]
  e$rec_index <- index[ie_str == "R"]
  e$icu_index <- index[ie_str == "ICU"]
  e$hosp_index <- index[ie_str == "H"]
  
  #### indices for calculating infections, hospitalizations, and deaths (code was originally in solve_model.R)
  e$v_ind_inf <- index[ie_str %in% c(v_exp_str, v_asym_inf_str, v_inf_str, "H", "ICU")]
  e$v_ind_hosp <- index[ie_str %in% c("H", "ICU")]
  e$v_ind_hosp_60p <- index[ie_str %in% c("H", "ICU") & ia %in% c(6:8)]
  e$v_ind_death <- index[ie_str %in% c("D")]
  
  #### matrices that have only 0s and 1s turn on and off school closure and social distancing among elderly in mixing matrix
  m_tmp <- matrix(0, ncol = ncol(mixing_matrix), nrow = nrow(mixing_matrix))
  ## school closure (assuming only affecting the population aged < 10)
  m_school <- m_tmp
  m_school[1, 1] <- 1
  e$m_school <- m_school
  
  #### age group mixing matrix
  ## Group index: 
  ## 1: 0-19
  ## 2: 20-39
  ## 3: 40-59
  ## 4: 60+
  tmp_ls <- list(c(1:2), c(3:4), c(5:6), c(7:9))
  m <- lapply(c(1:length(tmp_ls)), function(x) {
    out_ls <- lapply(c(1:length(tmp_ls)), function(y, z = x) {
      out_mat <- m_tmp 
      out_mat[tmp_ls[[z]], tmp_ls[[y]]] <- 1
      out_mat
    })
    names(out_ls) <- paste0(c(1:length(out_ls)))
    return(out_ls)
  })
  names(m) <- paste0("m", c(1:length(m)))
  e$m_60p <- m$m4[[4]]
  e$m_60p_w_young <- m$m4[[1]] + m$m4[[2]] + m$m4[[3]] + 
    m$m1[[4]] + m$m2[[4]] + m$m3[[4]]
  
  m <- unlist(m, recursive = F)
  names(m) <- gsub("[.]", "_", names(m))
  e$m <- m
  
  return(e)
}


#' Generate Indices
#'
#' takes vector of epigroups, comorbidity groups, and age groups and generates the indices that will reference the matching state in the
#' vector of model states. Note that the number of exposed and infectious compartments used in the model is a variable that is specified in
#' the file 'parameters_v3.R'. The get_ind function is called in the parameters() function to create a data.frame for indexing that is used throughout
#' the rest of the model.
#' @param ncg total number of comorbidity groups
#' @param neg total number of epigroups
#' @param nag total number of age groups
#' @param epi_groups a vector of epi states
#'
#' @return
#' a list of vectors of integers which reference the matching state in the vector of model states
#' 
#' @import plyr
#' 
#' @export
get_ind <- function(ncg = 2, neg = 61, nag = 9, epi_groups){
  
  df_ind <- expand.grid(ie = 1:neg, ic = 0:(ncg - 1), ia = 0:(nag - 1))
  df_ind$ie_str <- mapvalues(x = df_ind$ie,
                             from = c(1:neg),
                             to = epi_groups)
  
  eg <- df_ind[, "ie"] 
  cg <- df_ind[, "ic"]
  ag <- df_ind[, "ia"]
  
  if( any(eg < 1 | eg > neg) ) { stop('epigroup index out of bounds') }
  if( any(cg < 0 | cg > ncg) ) { stop('comorbidity index out of bounds') }
  if( any(ag < 0 | ag > nag) ) { stop('age index out of bounds') }
  
  index <- eg * ncg + cg + ag * neg * ncg - (ncg - 1)
  
  return(list(ie = df_ind$ie, ic = df_ind$ic, ia = df_ind$ia, 
              ie_str = df_ind$ie_str, index = index))
}
