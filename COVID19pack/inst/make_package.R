# Set working directory in `COVID19pack/` folder
# setwd("..")
main_dir  <- file.path(dirname(getwd()))
proj_path <- paste0(main_dir, "/COVID19proj")
data_path <- paste0(main_dir, '/data/')
pack_path <- paste0(main_dir, '/COVID19pack')

# Need to change the path!!!!!!!!!
sapply(c("mixing_matrix_MN_sym.rda", 
         "2020_03_23_MN_Init_Cases.rda", 
         "MN_population_data.rda"), 
       function(x) { load(paste0(data_path, x), envir = globalenv()) })

setwd(pack_path)
usethis::use_data(m_init_cases, overwrite = T)
usethis::use_data(mixing_matrix, overwrite = T)
usethis::use_data(population, overwrite = T)

devtools::document()
package_loc <- devtools::build()
install.packages(package_loc, repos = NULL) #, lib = proj_path)

#### Removing tar.gz file
# setwd("..")
# file.remove("COVID19pack_1.0.tar.gz")
