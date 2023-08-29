# This file is to generate results of sensitivity analysis.

rm(list = ls())

# Loading COVID-19 package
library(COVID19pack)

#### Set working directory
setwd(paste0(getwd(), "/COVID19proj"))
# This file includes scn_sim function, to calculate results for sensitivity analysis of 
# risk of transmissibility (beta), relative risk of hospitalizations and different level of 
# fatigue contact reductions (shown in Figure 3)
source("code/scn_fun.R") 
# This file includes scntri_sim function, to calculate results for sensitivity analysis of 
# triggering threshold, adjustment time and NPI effectiveness (shown in Figure S3)
source("code/scntri_sim.R")


# Basecase
scn_sim(name="BaseCase")

# Fatigue Contact Reductions
# -0.05
scn_sim(name = "CR05", fatigue = 0.05)
# -0.15
scn_sim(name = "CR15", fatigue = 0.15)
# -0.20
scn_sim(name = "CR20", fatigue = 0.20)
# -0.075
scn_sim(name = "CR075", fatigue = 0.075)
# -0.125
scn_sim(name = "CR125", fatigue = 0.125)
# -0.175
scn_sim(name = "CR175", fatigue = 0.175)

# Beta After 9/1 change
# x1.2
scn_sim(name = "betax12", beta_change = TRUE, delta_beta_after = 1.2)
# x1.5
scn_sim(name = "betax15", beta_change = TRUE, delta_beta_after = 1.5)
# x2.0
scn_sim(name = "betax2", beta_change = TRUE, delta_beta_after = 2.0)
# x0.8
scn_sim(name = "betax08", beta_change = TRUE, delta_beta_after = 0.8)

# Relative risk of hosp
# x0.5
scn_sim(name = "hosp05", hosp_change = TRUE, hosp_mul = 0.5)
# x0.7
scn_sim(name = "hosp07", hosp_change = TRUE, hosp_mul = 0.7)
# x0.9
scn_sim(name = "hosp09", hosp_change = TRUE, hosp_mul = 0.9)
# x1.2
scn_sim(name = "hosp12", hosp_change = TRUE, hosp_mul = 1.2)



# NPI Effectiveness
# Effectiveness 0.8 (baseline)
scntri_sim(name = "eff8")
# Effectiveness 0.9 (Scenario E)
scntri_sim(name = "eff9", eff = 0.9)
# Effectiveness 0.85 
scntri_sim(name = "eff85", eff = 0.85)
# Effectiveness 0.75 
scntri_sim(name = "eff75", eff = 0.75)
# Effectiveness 0.7 
scntri_sim(name = "eff7", eff = 0.7)
# Effectiveness 0.6 
scntri_sim(name = "eff6", eff = 0.6)
# Effectiveness 0.65
scntri_sim(name = "eff65", eff = 0.65)

## Triggering threshold
# Threshold 8 (baseline)
scntri_sim(name = "thr8")
# Threshold 6
scntri_sim(name = "thr6", thresh = 6)
# Threshold 4
scntri_sim(name = "thr4", thresh = 4)
# Threshold 10
scntri_sim(name = "thr10", thresh = 10)
# Threshold 12
scntri_sim(name = "thr12", thresh = 12)
# Threshold 6.5
scntri_sim(name = "thr65", thresh = 6.5)
# Threshold 7
scntri_sim(name = "thr7", thresh = 7)
# Threshold 7.5
scntri_sim(name = "thr75", thresh = 7.5)
# Threshold 9
scntri_sim(name = "thr9", thresh = 9)
# Threshold 8.5
scntri_sim(name = "thr85", thresh = 8.5)

## Adjustment time
# Adjustment time (baseline:7 (days))
scntri_sim(name = "adj7")
# Adjustment time 3 (days)
scntri_sim(name = "adj3",grace_period = 3)
# Adjustment time 5 (days)
scntri_sim(name = "adj5",grace_period = 5)
# Adjustment time 14 (days)
scntri_sim(name = "adj14",grace_period = 14)
# Adjustment time 21 (days)
scntri_sim(name = "adj21",grace_period = 21)

