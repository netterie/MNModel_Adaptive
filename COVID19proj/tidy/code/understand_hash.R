

# Try to understand the numeric indexes that were created
library(plyr)

# The inner workings of "gen_ind()", defined in "gen_states_index.R"
if (1==0) {
  comorbidity_groups = c("c0", "c1")
  ncg <- 2
  n_exposed_states <- 2
  n_infected_states <- 3
  
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
  age_groups <- seq(0,80,by=10)
  n_age_groups <- length(age_groups)
  
  ## Names for output objects
  comorbidity_groups <- c("c0", "c1")
  nam <- c(outer(outer(comorbidity_groups, epi_groups, FUN = "paste"), age_groups , FUN = "paste"))
  
  ## Create index list
  # ls_ix <- get_ind(ncg = ncg, neg = n_epi_groups, nag = n_age_groups, epi_groups)
  
  ncg = ncg
  neg = n_epi_groups
  nag = n_age_groups
  
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
  
  epi_states_index <- list()
  epi_states_index$ls_ix <- list(ie = df_ind$ie, ic = df_ind$ic, ia = df_ind$ia,
              ie_str = df_ind$ie_str, index = index)
  
  
  
}

if (1==0) {
  
  # Inner workings of create_hash_table 
  
  # Extract numbers of states, strata, infected and exposed states
  #nag <- epi_states_index$n_age_groups
  #ncg <- epi_states_index$ncg
  nes <- 2
  nis <- 3
  
  # From gen_states_index
  v_exp_str <- paste0("E", 1:nes)
  v_qe_str <- paste0("QE",1:nes)
  v_inf_str <- paste0("I", 1:nis)
  v_asym_inf_str <- paste0("AI", 1:nis)
  v_qi_str <- paste0("QI", 1:nis)
  v_qai_str <- paste0("QAI", 1:nis)
  epi_groups <- c("S", v_exp_str, v_asym_inf_str, v_inf_str,
                  "QS",v_qe_str, v_qai_str, v_qi_str, 
                  "H", "ICU", "R","RD","D", "HD", "CH")
  epi_groups_ls <- list(v_exp_str = v_exp_str,
                        v_qe_str = v_qe_str,
                        v_inf_str = v_inf_str,
                        v_asym_inf_str = v_asym_inf_str,
                        v_qi_str = v_qi_str,
                        v_qai_str = v_qai_str) 
  
  epi_states_index$epi_groups_ls <- epi_groups_ls
  
  ie <- epi_states_index$ls_ix$ie
  ic <- epi_states_index$ls_ix$ic
  ia <- epi_states_index$ls_ix$ia
  ie_str <- epi_states_index$ls_ix$ie_str
  index <- epi_states_index$ls_ix$index

  # Name the string ID vectors in epi_groups_ls:
  # $v_exp_str
  # [1] "E1" "E2"
  # $v_qe_str
  # [1] "QE1" "QE2"
  # $v_inf_str
  # [1] "I1" "I2" "I3"
  # $v_asym_inf_str
  # [1] "AI1" "AI2" "AI3"
  # $v_qi_str
  # [1] "QI1" "QI2" "QI3"
  # $v_qai_str
  # [1] "QAI1" "QAI2" "QAI3"
  for (i in c(1:length(epi_states_index$epi_groups_ls))) {
    assign(names(epi_states_index$epi_groups_ls)[i], epi_states_index$epi_groups_ls[[i]])
  }
    
  e <- new.env()
  
  #### indices for difference vector "d" in covid_19_model_function_v3
  
  #JKB Note: create a data frame mapping ie_str to index, 
  # to help decode some of these indexes 
  decode <- data.frame(id=ie_str, ix=index)
  
  #JKB Note: I am reorganizing these
  
  # Each type of state: S, E, AI, I, R, D, H, ICU, and Q versions
  # Plus RD, HD and CH
  e$v_S_ind <- index[ie_str == "S"]
  e$v_QS_ind <- index[ie_str == "QS"]
  e$v_E_ind <- index[ie_str %in% v_exp_str]
  e$v_AI_ind <- index[ie_str %in% v_asym_inf_str]
  e$v_I_ind <- index[ie_str %in% v_inf_str]
  e$v_QS_ind <- index[ie_str == "QS"]
  e$v_QE_ind <- index[ie_str %in% v_qe_str]
  e$v_QAI_ind <- index[ie_str %in% v_qai_str]
  e$v_QI_ind <- index[ie_str %in% v_qi_str]
  e$v_H_ind <- index[ie_str == "H"]
  e$v_ICU_ind <- index[ie_str == "ICU"]
  e$v_R_ind <- index[ie_str == "R"]
  e$v_RD_ind <- index[ie_str == "RD"]
  e$v_D_ind <- index[ie_str == "D"]
  e$v_HD_ind <- index[ie_str == "HD"]
  e$v_CH_ind <- index[ie_str == "CH"]
  
  # More states - comment said they were 
  # "for lambda calculation in lambda.R file"
  # All D states (This one is a duplicate)
  e$v_ind_D_by_age <- index[ie_str == "D"]
  # All I states (AI and I)
  e$v_ind_inf_by_age <- index[ie_str %in% c(v_asym_inf_str, v_inf_str)]
    # JKB: Code used to decode
    unique(ie_str[ie_str %in% c(v_asym_inf_str, v_inf_str)])
  # All QI states (AI and I)
  e$v_ind_quar_inf_by_age <- index[ie_str %in% c(v_qai_str, v_qi_str)]
    unique(ie_str[ie_str %in% c(v_qai_str, v_qi_str)])
  # All Q states (QS, QAI, QR)
  e$v_ind_quar_noninf_by_age <- index[ie_str %in% c("QS",v_qe_str,"QR")]
  
  # More states, ALL DUPLICATES - comment said they were 
  # "indices for increase in testing (originally in covid_19_model_function())"
  e$s_index <- index[ie_str == "S"]
  # All E states
  e$exp_index <- index[ie_str %in% v_exp_str] 
    unique(ie_str[ie_str %in% v_exp_str])
  # All I states
  e$inf_index <- index[ie_str %in% v_inf_str]
    unique(ie_str[ie_str %in% v_inf_str])
  # All AI statues
  e$asym_index <- index[ie_str %in% v_asym_inf_str]
    unique(ie_str[ie_str %in% v_asym_inf_str])
  e$rec_index <- index[ie_str == "R"]
  e$icu_index <- index[ie_str == "ICU"]
  e$hosp_index <- index[ie_str == "H"]
  
  
  # "E1" states - i.e., the first E state
  e$v_E_ind1 <- e$v_E_ind[seq(1, length(e$v_E_ind), nes)]
    # JKB: Code used to decode
    tmp <- df_ix[decode$ix %in% e$v_E_ind,]
    tmp[seq(1, length(e$v_E_ind), nes),]
  
  # This is also "E1" - i.e., all "E" states except the last one?
  e$v_E_ind_rm_last <- e$v_E_ind[setdiff(c(1:length(e$v_E_ind)), seq(nes, length(e$v_E_ind), nes))]
    tmp[setdiff(c(1:length(e$v_E_ind)), seq(nes, length(e$v_E_ind), nes)),]
  
  # "E2" states - i.e., all "E" states except the first one?
  e$v_E_ind_rm_first <- e$v_E_ind[setdiff(c(1:length(e$v_E_ind)), seq(1, length(e$v_E_ind), nes))]
    tmp[setdiff(c(1:length(e$v_E_ind)), seq(1, length(e$v_E_ind), nes)),]
  
  # "E2" states - i.e, the last "E" state
  e$v_E_ind_last <- e$v_E_ind[seq(nes, length(e$v_E_ind), nes)]
    cbind(tmp[seq(nes, length(e$v_E_ind), nes),], tmp[tmp$id=="E2",])
  
  # Repeat the process for AI, I, QE, QAI, and QI: first, all but last, all but first, and last
  e$v_AI_ind1 <- e$v_AI_ind[seq(1, length(e$v_AI_ind), nis)]
  e$v_AI_ind_rm_last <- e$v_AI_ind[setdiff(c(1:length(e$v_AI_ind)), seq(nis, length(e$v_AI_ind), nis))]
  e$v_AI_ind_rm_first <- e$v_AI_ind[setdiff(c(1:length(e$v_AI_ind)), seq(1, length(e$v_AI_ind), nis))]
  e$v_AI_ind_last <- e$v_AI_ind[seq(nis, length(e$v_AI_ind), nis)]
  
  e$v_I_ind1 <- e$v_I_ind[seq(1, length(e$v_I_ind), nis)]
  e$v_I_ind_rm_last <- e$v_I_ind[setdiff(c(1:length(e$v_I_ind)), seq(nis, length(e$v_I_ind), nis))]
  e$v_I_ind_rm_first <- e$v_I_ind[setdiff(c(1:length(e$v_I_ind)), seq(1, length(e$v_I_ind), nis))]
  e$v_I_ind_last <- e$v_I_ind[seq(nis, length(e$v_I_ind), nis)]
  
  e$v_QE_ind1 <- e$v_QE_ind[seq(1, length(e$v_QE_ind), nes)]
  e$v_QE_ind_rm_last <- e$v_QE_ind[setdiff(c(1:length(e$v_QE_ind)), seq(nes, length(e$v_QE_ind), nes))]
  e$v_QE_ind_rm_first <- e$v_QE_ind[setdiff(c(1:length(e$v_QE_ind)), seq(1, length(e$v_QE_ind), nes))]
  e$v_QE_ind_last <- e$v_QE_ind[seq(nes, length(e$v_QE_ind), nes)]
  
  e$v_QAI_ind1 <- e$v_QAI_ind[seq(1, length(e$v_QAI_ind), nis)]
  e$v_QAI_ind_rm_last <- e$v_QAI_ind[setdiff(c(1:length(e$v_QAI_ind)), seq(nis, length(e$v_QAI_ind), nis))]
  e$v_QAI_ind_rm_first <- e$v_QAI_ind[setdiff(c(1:length(e$v_QAI_ind)), seq(1, length(e$v_QAI_ind), nis))]
  e$v_QAI_ind_last <- e$v_QAI_ind[seq(nis, length(e$v_QAI_ind), nis)]
  
  e$v_QI_ind1 <- e$v_QI_ind[seq(1, length(e$v_QI_ind), nis)]
  e$v_QI_ind_rm_last <- e$v_QI_ind[setdiff(c(1:length(e$v_QI_ind)), seq(nis, length(e$v_QI_ind), nis))]
  e$v_QI_ind_rm_first <- e$v_QI_ind[setdiff(c(1:length(e$v_QI_ind)), seq(1, length(e$v_QI_ind), nis))]
  e$v_QI_ind_last <-  e$v_QI_ind[seq(nis, length(e$v_QI_ind), nis)]
  
  
  #### indices for calculating infections, hospitalizations, and deaths (code was originally in solve_model.R)
  # All E, AI, I, H and ICU states
  e$v_ind_inf <- index[ie_str %in% c(v_exp_str, v_asym_inf_str, v_inf_str, "H", "ICU")]
    unique(ie_str[ie_str %in% c(v_exp_str, v_asym_inf_str, v_inf_str, "H", "ICU")])
  # H and ICU states
  e$v_ind_hosp <- index[ie_str %in% c("H", "ICU")]
  # I and ICU states for 60-plus
  e$v_ind_hosp_60p <- index[ie_str %in% c("H", "ICU") & ia %in% c(6:8)]
  # Death - again, this is a duplicate
  e$v_ind_death <- index[ie_str %in% c("D")]
  
  #### matrices that have only 0s and 1s turn on and off school closure and social distancing among elderly in mixing matrix
  # JKB: do a little setup to get mixing_matrix loaded 
  # JKB: this should be loaded as an input matrix, not created here. It's a matrix of zeros for every 
  # age-age interaction except for a 1 entry in [1,1] (ages 0-10 interacting with ages 0-10)
  # Not sure why school closure would only affect those under aged 10
  source(file.path("COVID19proj", "tests", "test_scenarios.R"))
  parms <- test_scenarios(nsamp=1,runmodel=FALSE)
  mixing_matrix <- parms$mixing_matrix
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
  
  # JKB: "m" is a list of lists. Each list is an age-age contact matrix of zeros
  # except for a small block of 1's for each of the 4 groups indicated above:
  # m$m1 - Group 1 with Groups 1,2,3,4
  # m$m2 - Group 2 with Groups 1,2,3,4
  # m$m3 - Group 3 with Groups 1,2,3,4
  # m$m4 - Group 4 with Groups 1,2,3,4
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
