# Set working directory in `COVID19pack/` folder
# setwd("..")
lib_path <- paste0(getwd(), "/COVID19proj")
setwd("COVID19pack/")

# Need to change the path!!!!!!!!!
tmp_path <- "Your Path/MNCOVID19_Trigger/data/"
sapply(c("mixing_matrix_MN_sym.rda", 
         "2020_03_23_MN_Init_Cases.rda", 
         "MN_population_data.rda"), 
       function(x) { load(paste0(tmp_path, x), envir = globalenv()) })

usethis::use_data(m_init_cases, overwrite = T)
usethis::use_data(mixing_matrix, overwrite = T)
usethis::use_data(population, overwrite = T)

devtools::document()
package_loc <- devtools::build()
install.packages(package_loc, repos = NULL) #, lib = lib_path)

#### Removing tar.gz file
# setwd("..")
# file.remove("COVID19pack_1.0.tar.gz")
