# COVID19pack

You could use the functions created for COVID19pack by sourcing the functions in the `R/` folder or build the package.

## Build COVID19pack

-   Run the `make_package.R` script in the `inst/` folder to generate the package. To create the package, the working directory has to be under `COVID19pack/`.

-   In addition, this package has three internal datasets `m_init_cases.rda`, `mixing_matrix.rda`, and `population.rda`. To create these datasets, remember to point to the proper data folder where stores the original datasets. The data folder that contains the original datasets should be under `MNCOVID19/` directory. Modify the `tmp_path` in line 7 to the proper path.
