
# ****************************************************************************
# Code to understand the initialize population code in the 
# get_initial_conditions() function, which is located in COVID19pack/R 
# and called from parameters.R. The two inputs are 
# 1. "parms" and 
# 2. "m_init_cases", which is in the COVID19pack/data folder
# ****************************************************************************

# Setup
library(tidyverse)
library(openxlsx)
this_folder <- file.path("COVID19proj", "tidy", "inputs")

# Read in parameter descriptions (copy-pasted from the help file)
parms_descriptions <- read.csv(file.path(this_folder, "parms_Adaptive_descriptions.csv"))

# Set up parms 
source(file.path("COVID19proj", "tests", "test_scenarios.R"))
parms <- test_scenarios(nsamp=1,runmodel=FALSE)

# ----------------------------------------------------------------------------
# m_init_cases
# ----------------------------------------------------------------------------

# This is a df-turned-matrix with one row per the following states: I, H, ICU, D and R.
# Columns are age groups, titled "age_grp1"...."age_grp9"
  m_init_cases <- as.matrix(m_init_cases)


# ----------------------------------------------------------------------------
# Step into get_initial_conditions
# ----------------------------------------------------------------------------

# get_initial_conditions <- function(parms, m_init_cases){

# ----------------------------------------------------------------------------
# Initialize some objects from parms
# ----------------------------------------------------------------------------

# Scalar that has a default of 0.019 set in parameters.R, but for some reason is .0141 in 
# parms_Adaptive_guide.xlsx, and is now 0.0297943
  init_cases_detected <- parms$init_cases_detected 

# MN population
  N <- parms$N
  
# Proportion in each age group
  age_prop <- parms$age_prop
# JKB: one entry per age group; defaults sum to 1.777. Within each age group, percent having a comorbidity?
# So, conditional on age, percent in the strata?
  comorbidity_prop_by_age <- parms$comorbidity_prop_by_age
# Keys for states E1/E2, I1/I2/I3, and AI1/AI2/AI3
  v_exp_str <- parms$epi_groups_ls$v_exp_str
  v_inf_str <- parms$epi_groups_ls$v_inf_str
  v_asym_inf_str <- parms$epi_groups_ls$v_asym_inf_str 
# Conditional on age, proportion asymptomatic
  prop_asym <- parms$prop_asymptomatic
# JKB: Conditional on age, proportion infected?
  prop_inf_by_age <- parms$prop_inf_by_age
# Number of infected and exposed state
  nis <- parms$n_infected_states
  nes <- parms$n_exposed_states
  
# ----------------------------------------------------------------------------
# Initialize compartment indexes from parms
# ----------------------------------------------------------------------------
  
# Hash indexes: vectors of length 450
  # Comorbidity group (0/1)
  ic <- parms$ls_ix$ic
  # Age group (0,1,...8)
  ia <- parms$ls_ix$ia
  # State ID: S, E1, E2, A1..A3, S1..S3, Q-states for all the aformentioned,
  # H, ICU, R, RD, D, HD, CH
  ie_str <- parms$ls_ix$ie_str
  # Numeric index 1...450 but NOT in numeric order due to the formula
  # used to construct it in set_states_eo.R. It counts up by 50s using the odds a
  # and then the evens, e.g. 1,3,5...49, 2,4,6....50, 51,53,55.....99, 52,54...100...
  # all the way up to 450. This sets the odds as comorbidity group==0 and the 
  # evens as comorbidity group==1. 
  index <- parms$ls_ix$index
  
# ----------------------------------------------------------------------------
# Estimate undetected cases
# ----------------------------------------------------------------------------
  # estimate undetected cases: scalar = 7587.305
  # Derived by summing all I, H and ICU cases per age group, then multiplying that
  # by proportion of cases undetected to get undetected cases by age. Sum over age 
  # and divide by proportion undetected to get total cases undetected
  # This could be more readable if interim calculations were assigned a good object
  # name, e.g. v_n_cases, n_cases_detected, etc.
  total_v_init_NDinf <- sum(colSums(m_init_cases) * (1 - init_cases_detected) / init_cases_detected)
  
  # Undetected cases by age 
  v_init_NDinf <- total_v_init_NDinf * prop_inf_by_age
  
# ----------------------------------------------------------------------------
# Initialize some vectors as zero
# ----------------------------------------------------------------------------
  v_init_NDinf_asym <- rep(0, 9)
  v_init_NDinf_sym <- rep(0, 9)
  
  # Initial starting vector
  # A vector with one slot per state-compartment (length 450)
  init_vec <- 0 * vector(mode = "numeric", length = length(index))
  #print(length(init_vec))
  
# ----------------------------------------------------------------------------
# Total cases by age (detected + undetected)
# ----------------------------------------------------------------------------
  # First, distribute people across age and comorbidity compartments
  # indexing will be ia = (i-1) and cg=0 or 1
  # cases (known and unknown by age)
  
  # This is a vector with the count of total cases (detected plus undetected) by age
  # Use of apply is confusing because the summation of cases across I,H, and ICU 
  # was already done using colSums when creating total_v_init_NDinf
  n_cases_by_age <- apply(m_init_cases, 2, sum) + v_init_NDinf

# ----------------------------------------------------------------------------
# Initialize state-compartments: S
# ----------------------------------------------------------------------------
  # Populate S state with everyone who is not a case, by strata
  ## Susceptible (eg=1) ##
  # Num in age group i with no co-morbidities
  init_vec[index[ie_str == "S" & ic == 0]] <- (N * age_prop - n_cases_by_age) * (1 - comorbidity_prop_by_age)
    
  # Num in age group i with at least one co-morbidity
  init_vec[index[ie_str == "S" & ic == 1]] <- (N * age_prop - n_cases_by_age) * comorbidity_prop_by_age
    
# ----------------------------------------------------------------------------
# Initialize state-compartments: I
# ----------------------------------------------------------------------------
  # Populate final I state with everyone who is a detected case, by strata
  
  ## Detected cases ##
  # Infected -> go to last infected state
  # NOTE: m_init_cases must be a matrix for this to work correctly
  # if it is a data.frame, init_vec will turn into a list
  # into a list? I'm going to use [[]] to keep the original structure
  # Hash needed: last I state, reduced by comorbidity status
  # I think this can be generalized using a loop, or requiring this
  # to be data analysis done ahead of here
  init_vec[index[ie_str == v_inf_str[length(v_inf_str)] & ic == 0]] <-
    m_init_cases["I", ] * (1 - comorbidity_prop_by_age)
  
  init_vec[index[ie_str == v_inf_str[length(v_inf_str)] & ic == 1]] <-
    m_init_cases["I", ] * comorbidity_prop_by_age 
  
# ----------------------------------------------------------------------------
# Initialize state-compartments: H
# ----------------------------------------------------------------------------
  # Populate H state with everyone who is a detected case, by strata
    
  # Hospitalized -> go to H (eg=16)
  init_vec[index[ie_str == "H" & ic == 0]] <- m_init_cases["H", ] * (1 - comorbidity_prop_by_age)
  init_vec[index[ie_str == "H" & ic == 1]] <- m_init_cases["H", ] * comorbidity_prop_by_age
  
# ----------------------------------------------------------------------------
# Initialize state-compartments: ICU
# ----------------------------------------------------------------------------
  # Populate ICU state with everyone who is a detected case, by strata
    
  # ICU -> go to ICU (eg=17)
  init_vec[index[ie_str == "ICU" & ic == 0]] <- m_init_cases["ICU", ] * (1 - comorbidity_prop_by_age)
  init_vec[index[ie_str == "ICU" & ic == 1]] <- m_init_cases["ICU", ] * comorbidity_prop_by_age
    
# ----------------------------------------------------------------------------
# Initialize state-compartments: R, D
# ----------------------------------------------------------------------------
  # Populate R state with everyone who is a detected case, by strata. Wouldn't this always be zero?
  
  # Recovered -> go to R (eg=18)
  init_vec[index[ie_str == "R" & ic == 0]] <- m_init_cases["R", ] * (1 - comorbidity_prop_by_age)
  init_vec[index[ie_str == "R" & ic == 1]] <- m_init_cases["R", ] * comorbidity_prop_by_age
  
  # Populate R state with everyone who is a detected case, by strata. Wouldn't this always be zero?
  
  # Dead --> got to D (eg=19)
  init_vec[index[ie_str == "D" & ic == 0]] <- m_init_cases["D", ] * (1 - comorbidity_prop_by_age)
  init_vec[index[ie_str == "D" & ic == 1]] <- m_init_cases["D", ] * comorbidity_prop_by_age

## Distribute undetected cases proportionately across exposed, asymptomatic, and infectious compartments
# ----------------------------------------------------------------------------
# Initialize state-compartments: E, AI and SI (use undetected cases estimate)
# STEP 1: Determine # of undetected cases to assign to each compartment, 
# across all strata
# ----------------------------------------------------------------------------
  
  # Number of E+I and AI compartments
  # NOTE: Suggest using variables consistently, e.g. nis (# infected states) and nes (# exposed states) 
  # are used further down but not here. 
  n_comp_s_e <- length(c(v_exp_str, v_inf_str)) # n_comp_s_e = nes + nis = 5
  n_comp_a <- length(v_asym_inf_str)            # n_comp_a = 3
  
  # JKB: I'm not sure why prop_asym is not the right dimensions. In parameters.R it looks like it 
  # should be. Rather than beat my head over it, I'm just going to re-define it here.
  # It's a (# strata) x (# age group) matrix of p(asymptomatic | case, age)
  ncg=2
  n_age_groups=9
  prop_asymp_19 = 0.408
  prop_asym <- matrix(0, nrow = ncg, ncol = n_age_groups)
  tmp <- c(1, 1.113, 1.028, 0.944, 0.845, 0.718, 0.521, 0.437, 0.437) * prop_asymp_19  # Davies et al.
  prop_asym[1, ] <- prop_asym[2, ] <- tmp
    
  # Undetected cases = proportion of undetected cases * probability that they are asymptomatic/symptomatic
  # Multiply that by (# infected + # exposed = 5) states and then divide by 
  # ( 2*(# exposed states) + (# infected) = 7 ) .... ?
  # NOTE: Suggest initializing vectors and defining them at the top, or as the code arises
  # These were defined earlier but simply set to zero. Unnecessary?
  v_init_NDinf_asym <- v_init_NDinf * prop_asym[1, ] * (nis + nes) / (2 * nes + nis) 
  v_init_NDinf_sym <- v_init_NDinf * (1 - prop_asym[1, ] * (nis + nes) / (2 * nes + nis))
    
# ----------------------------------------------------------------------------
# Initialize state-compartments: E, AI and SI (use undetected cases estimate)
# STEP 2: Split undetected cases in each state across strata
# ----------------------------------------------------------------------------
  
  # Suggest organizing by state, not by strata
  
  # Question: can we delete the second line, adding the prior counts? Aren't those 
  # zero because we're initializing the population?
  
  # Question: is there a more systematic way to assign things by strata -- use
  # a loop, if things are repeatedly getting multiplied by the same strata vector?
  # Maybe a pre-step should be to get all strata proportions reduced to one master-strata vector,
  # ordered correctly according to index, and then go through and implement the strata-
  # related calculations using a loop
  
  # no comorbidities
  #asymptomatic
  # Divide the asymptomatic cases evenly across asymptomatic compartments,
  # then multiply by comorbidity proportion and repeat that value for each
  # of the state-compartments
  # NOTE: very hard-coded, needs generalizing
  init_vec[index[ie_str %in% c(v_asym_inf_str) & ic == 0]] <-
    init_vec[index[ie_str %in% c(v_asym_inf_str) & ic == 0]] +
    rep((v_init_NDinf_asym / (n_comp_a) * (1 - comorbidity_prop_by_age)), 
        each = length(v_asym_inf_str))
  
  # Divide the symptomatic cases evenly across exposed and symptomatic compartments,
  # then multiply by comorbidity proportion and repeat that value for each of the 
  # state-compartments
  # NOTE: very hard-coded, needs generalizing
  # Hash needed: E-SI
  #symptomatic and exposed
  init_vec[index[ie_str %in% c(v_exp_str, v_inf_str) & ic == 0]] <-
    init_vec[index[ie_str %in% c(v_exp_str, v_inf_str) & ic == 0]] +
    rep((v_init_NDinf_sym / (n_comp_s_e) * (1 - comorbidity_prop_by_age)), 
        each = length(c(v_exp_str, v_inf_str)))

  # comorbidities
  #asymptomatic
  init_vec[index[ie_str %in% c(v_asym_inf_str) & ic == 1]] <-
    init_vec[index[ie_str %in% c(v_asym_inf_str) & ic == 1]] +
    rep((v_init_NDinf_asym / (n_comp_a) * (comorbidity_prop_by_age)), 
        each = length(v_asym_inf_str))
  
  #symptomatic and exposed
  init_vec[index[ie_str %in% c(v_exp_str, v_inf_str) & ic == 0]] <-
    init_vec[index[ie_str %in% c(v_exp_str, v_inf_str) & ic == 0]] +
    rep((v_init_NDinf_sym / (n_comp_s_e) * (comorbidity_prop_by_age)), 
        each = length(c(v_exp_str, v_inf_str)))
  