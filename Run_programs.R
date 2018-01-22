
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
### To run analytical solution solutions for various model assumptions.
### All cases are labelled in numeric orders.
### The detailed descriptions of each case is available in their respective script comments
source("R/Run_analytical_solutions.R")






##### ------------------------ Clear workspace ---------------------------- #####
rm(list=ls(all=TRUE))

#### End


