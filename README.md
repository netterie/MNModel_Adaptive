# Accompanying R code for the forthcoming publication 
# Sanstead et al., Adaptive COVID-19 Mitigation Strategies: Trade-Offs Between Trigger Thresholds, Response Timing, and Effectiveness, MDM Policy & Practice, 2023.
## How to use:

First, go to COVID19pack/inst/make_package.R to install the COVID19 model package on your machine.

**JKB Note**: Open a new R session using the COVID19pack.Rproj I added within the CODIV19pack folder, and run make_package.R from there.
To run the analysis scripts below, exit that RProj, open the COVID19_Trigger.RProj, and work from there.

The following scripts may be of interest:

-   The primary script that is used to run the base case scenarios A-F is Scenarios.R
-   To re-run model calibration, use the script calibABC.R
-   To run sensitivity analyses, use Sensitivity_Analysis.R

Some scripts may take a long time to run. Take care to set your working directory appropriately.
