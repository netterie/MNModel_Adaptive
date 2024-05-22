# ****************************************************************************
# Code to experiment with tidying parameter inputs and prep
# ****************************************************************************

# Setup
library(tidyverse)
library(openxlsx)
this_folder <- file.path("COVID19proj", "tidy", "inputs")

# Read in parameter descriptions (copy-pasted from the help file)
parms_descriptions <- read.csv(file.path(this_folder, "parms_Adaptive_descriptions.csv"))

# Set up parms using the historic approach
source(file.path("COVID19proj", "tests", "test_scenarios.R"))
parms <- test_scenarios(nsamp=1,runmodel=FALSE)

# Save parms element names and element types to csv
parms_type <- sapply(parms, function(x) { class(x)[1]})
parms_length <- sapply(parms, function(x) { length(x) })
parms_value <- sapply(parms, function(x) {
  if (length(x)==1 & !is.list(x)) {
    return(as.numeric(x))
  } else return(NA)
})
return_dim <- function(x, d) {
  thisdim = dim(x)[d]
  if (!is.null(thisdim)) {
    return(thisdim) 
  } else if (is.list(x)) {
    return(length(x[[1]])[1])
  } else return(NA)
}
parms_dim1 <- sapply(parms, return_dim, 1) 
parms_dim2 <- sapply(parms, return_dim, 2)
parms_info <- data.frame(Name=names(parms_type), 
                         Type=parms_type,
                         Length=parms_length,
                         Value=parms_value,
                         Dim1=parms_dim1,
                         Dim2=parms_dim2)

# Merge in descriptions
parms_guide <- left_join(parms_info, parms_descriptions, by="Name")

# Save
write.csv(parms_guide, file.path(this_folder, "parms_Adaptive_guide.csv"), row.names=FALSE, na="")

# Save matrices
parms_mats_names <- names(parms_type[parms_type=="matrix"])
parms_mats <- parms[parms_mats_names]
parms_mats_shortnames <- gsub("reduction", "red", parms_mats_names)
names(parms_mats) <- parms_mats_shortnames
write.xlsx(parms_mats, file.path(this_folder, "parms_Adaptive_matrices.xlsx"))

# Save vectors - extract from lists, then add single vectors
# Vectors from lists
parms_vector_names <- c(
  "prop_hosp_ls",
  "prop_ICU_ls",
  "p_h_exit",
  "p_icu_exit",
  "prob_hosp_death",
  "prob_icu_death"
)
parms_vectors <- lapply(parms_vector_names, 
                        function(x) {
                          cat(x, "...")
                          # Get the list object
                          y = parms[[x]]
                          # For each list element, return a matrix with names
                          ret <- sapply(y, function(z) z)
                          ynames <- names(y)
                          if (is.null(ynames)) ynames <- 1:ncol(ret)
                          colnames(ret) <- paste(x, ynames, sep="_")
                          return(ret)
                        })

# Add single vectors
parms_vector_names <- c(
  "prop_hosp",
  "prop_ICU", 
  "prob_icu_death_no_bed",
  "prop_inf_die",
  "prop_asymptomatic"
)
parms_vectors2 <- sapply(parms_vector_names, function(x) parms[[x]])
parms_vectors_df <- cbind(do.call(cbind, parms_vectors), parms_vectors2)


write.csv(parms_vectors_df, file.path(this_folder, "parms_Adaptive_vectors.csv"), 
          row.names=FALSE)
