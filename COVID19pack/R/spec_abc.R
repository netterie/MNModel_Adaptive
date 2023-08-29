# ### specification of calibration parameters 
# 

#' The start_date is always set on 2020-03-22
#' @export
set_start_date <- function() {
  as.Date("2020-03-22")
}

#' @keywords internal
spec_parm23 <- function() {
  v_names_params <- c("init_cases_detected",
                      "p_h_50_59",
                      "p_h_60p",
                      "p_dying_home_80",
                      "prop_asymp_19",
                      "beta",
                      "weight_60p_init",
                      "social_distancing_contact_reduction",
                      paste0("sip_cr", c(1:5)),
                      "prop_of_sip1", 
                      "beh_cr1", 
                      "prop_of_sip2", 
                      "beh_cr2")
  v_lb <- c(0.001, 0, 0, 0.000000001, 0.20, 0.001, 0.01, 0.1,
            rep(0, 5), rep(0, 4))
  v_ub <- c(0.25, 0.4, 0.4, 0.15, 0.6, 0.04, 0.9, 0.5,
            rep(1, 5), rep(1, 4))
  
  parm_bd_df <- data.frame(parm_names = v_names_params,
                           lb = v_lb,
                           ub = v_ub)
  return(parm_bd_df)
}

#' @keywords internal
gen_abc_param23 <- function(x, epi_states_index, hash_table,
                            mask_cr_mat = NULL, # estimated contact reduction matrix including masking
                            start_time_mask_policy = Inf,
                            beta_after = 0.02, 
                            beta_change = FALSE,
                            hosp_change = FALSE,
                            hosp_mul = 1) {
  
  sip_cr_mat <- matrix(c(x[[9]],  x[[10]], x[[10]], x[[11]],
                         x[[10]], x[[12]], x[[12]], x[[11]],
                         x[[10]], x[[12]], x[[12]], x[[11]],
                         x[[11]], x[[11]], x[[11]], x[[13]]),
                       nrow = 4, byrow = T)
  
  beh_cr_mat <- matrix(c(x[[14]]*x[[9]], x[[15]],         x[[15]],         x[[17]],
                         x[[15]],        x[[16]]*x[[12]], x[[16]]*x[[12]], x[[17]],
                         x[[15]],        x[[16]]*x[[12]], x[[16]]*x[[12]], x[[17]],
                         x[[17]],        x[[17]],         x[[17]],         x[[13]]),
                       nrow = 4, byrow = T)
  
  parameters(epi_states_index = epi_states_index,
             hash_table = hash_table,
             init_cases_detected = x[[1]],
             p_h_50_59 = c(x[[2]], x[[2]]),
             p_h_60p = c(x[[3]], x[[3]]),
             p_dying_home_80 = x[[4]],
             prop_asymp_19 = x[[5]],
             beta = x[[6]],
             weight_60p_init = x[[7]],
             social_distancing_contact_reduction = x[[8]],
             sip_cr_mat = sip_cr_mat,
             beh_cr_mat = beh_cr_mat,
             mask_cr_mat = mask_cr_mat, 
             start_time_mask_policy = start_time_mask_policy,
             beta_after = beta_after,
             beta_change = beta_change,
             hosp_mul = hosp_mul,
             hosp_change = hosp_change)
}
