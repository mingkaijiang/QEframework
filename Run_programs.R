
####################### Master Script for QE framework analysis ###############################
###
### Summary:
#### To generate analytical solutions based on quasi-equilibrium framework 
###  to test effect of alternative model assumptions on model behaviors
###
### Author: Mingkai Jiang (m.jiang@westernsydney.edu.au)
### 


#### ------------------------ General system stuffs ------------------------ #####
### Make sure everything is clear
rm(list=ls(all=TRUE))

### Get current date
date<-Sys.Date()

### pre-processing script, read in all R packages
source("QE_Analyses/prepare_R.R")

#### ------------------------ Run analytical stuffs ------------------------ #####
### To run analytical solution solutions for various model assumptions
### Currently 10 cases, for details see the comments in the source code
### Two dataframe are generated: equilDF and constraintDF
source("R/Run_analytical_solutions.R")

#### ------------- Checking effect of wood stoichiometry in analytical solution ----------- #####
source("R/Effect_of_wood_stoichiometry.R")


#### ------------- Plotting CO2 fertilization effect through time ----------- #####
source("R/CO2_fertilization_summary.R")


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


