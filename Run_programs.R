
#### Main program
###
### Summary:
### 1. To prepare GDAY met forcing data and parameter files
### 2. To run GDAY and output each simulations into the corresponding folders
### 3. To quality check GDAY simulation results and post-processing them
### 4. To generate analytical solutions for GDAY simulations
### 5. To prepare manuscript figures and tables 
###
### Author: Mingkai Jiang (m.jiang@westernsydney.edu.au)
### 
### Warning: large hard disk space needed,
### because the simulation creates many large files
### (but many simulation files will be deleted after processing)


#### ------------------------ General system stuffs ------------------------ #####
### Make sure everything is clear
rm(list=ls(all=TRUE))

### Get current date
date<-Sys.Date()

### read in all R packages
source("R/prepare_R.R")

#### ------------------------ Prepare GDAY stuffs -------------------------- #####
### Create met data for gday simulations
#source("GDAY/pre_processing/create_monthly_met_for_GDAY.R")
source("GDAY/pre_processing/create_daily_met_for_GDAY.R")

### compile gday program and send to simulation folders
source("GDAY/pre_processing/Make_GDAY_and_Send_To_Folders.R")

### Here need a script to modify the R scripts parameters for each simulations
source("GDAY/pre_processing/Paste_R_script_to_folders.R")


#### ------------------------ Run GDAY simulations ------------------------- #####
### Run GDAY simulations, using either the python or R wrapper file
### Current setting use R
### Only needs GDAY simulation to verify the analytical solution
### Currently only implemented GDAY simulations for Run 1 - 3
### i.e. CNP with variable wood, CN with variable wood, CNP with fixed wood stoichiometry
source("GDAY/pre_processing/Run_GDAY.R")


#### ------------------ Post-processing GDAY simulations ------------------- #####
### Convert from monthly to annual data and save to analyses subfolders
### This step is the only "must-run" step for post-processing purpose
#source("GDAY/post_processing/Convert_GDAY_monthly_to_annual.R")   
source("GDAY/post_processing/Convert_GDAY_daily_to_annual.R")   


### delete all raw GDAY output files as they are large
source("GDAY/post_processing/delete_raw_gday_files.R")

### Mass balance checks, takes very long to run!
# source("GDAY/post_processing/mass_balance.R")

### Plot time series spin up files for each simulations
source("GDAY/post_processing/Transient_spin_up_plot.R")

### Check spin-up and transient continuity
### co2_amb and co2_ele only pools, starting from last 10 years of equilibration
source("GDAY/post_processing/Check_continuity.R")

### plot continuity starting from transient year 1 and specify endyear for pools and flxues
### Only for elevated CO2 runs
source("GDAY/post_processing/Check_continuity_transient.R")

### plot gday simulated quasi-equil points under aCO2 and eCO2
### Note: need to specify years when L and VL equilibrates
###       better to consider an automatic process to pick these years
source("GDAY/post_processing/Plot_GDAY_quasi_equil_constraints.R")


#### ------------------------ Run analytical stuffs ------------------------ #####
### To run analytical solution solutions for various model assumptions
### Currently 10 cases, for details see the comments in the source code
### Two dataframe are generated: equilDF and constraintDF
source("R/Run_analytical_solutions.R")


#### ------------- Checking GDAY matches with analytical results ----------- #####
### generate a table for comparison of the equilibrium points
source("R/Check_analytical_gday_matches.R")

# Note: VL does not 100% matched with GDAY simulated VL, possible reasons:
# GDAY equilibrium definition: < 1e-2 difference might not mean equilibrium.
# Note: when NC and PC matched, NPP does not match, possible reason:
# the full photosynthesis model is simplified in the analytical equations (i.e. hourly vs. daily)

#### ------------- Checking effect of P limitation in analytical solution ----------- #####
source("R/Effect_of_P_limitation.R")

#### ------------- Checking effect of wood stoichiometry in analytical solution ----------- #####
source("R/Effect_of_wood_stoichiometry.R")

#### ------------- Checking effect of P limitation in GDAY key variables ----------- #####
source("R/Effect_of_P_GDAY.R")

#### ------------- Plotting CO2 fertilization effect through time ----------- #####
source("R/CO2_fertilization_summary.R")

#### ------------- Case studies - EucFACE and AmazonFACE ----------- #####
#source("R/Effect_of_P_limitation_EucFACE.R")
#
#source("R/Effect_of_P_limitation_AmazonFACE.R")
#
#source("R/Run_analytical_solutions_EucFACE.R")
#
#source("R/CO2_fertilization_summary_EucFACE.R")
#
#source("R/Run_analytical_solutions_AmazonFACE.R")
#
#source("R/CO2_fertilization_summary_AmazonFACE.R")
#
#### ---------------- Generate manuscript figures and tables --------------- #####
### To generate manuscript figures
source("Plots/Figure_generating.R")

### To generate manuscript tables (or statistics used for generating the tables)
source("Tables/Table_generating.R")

### To generate a conceptual animated figure to show how things are moving within 
### the quasi-equilibrium graph
source("Plots/Animated_Figure_Generation.R")


##### ------------------------ Clear workspace ---------------------------- #####
rm(list=ls(all=TRUE))

#### End


