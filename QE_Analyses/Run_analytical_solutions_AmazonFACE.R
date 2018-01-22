
#### Run analytical functions for each sub-simulations of GDAY
#### Output the necessary dataframe onto a cohesive framework,
#### so that they are available for cross-match with GDAY simulation results
####
#### This is the main script to call individual analytical run scripts and functions
####
#### Run definitions:
#### Run 1: baseline, variable wood stoichiometry, N and P cycles on,
####        implicit mineral N and P pools
#### Run 2: same as Run1, but P cycle off

################################################################################
#### Create dataframes to store all the data
constraintDF_FACE <- create_FACE_constraint_DF()
equilDF_FACE <- create_FACE_equil_DF()

#### Step 1: simply run analytical solution and plot quasi-equil plots
### f.flag: = 1 simply plot analytical solution graph
###         = 2 return constraintDF
###         = 3 return equilDF
Perform_Analytical_Run1_AmazonFACE(f.flag = 1, constraintDF_FACE, equilDF_FACE)
Perform_Analytical_Run2_AmazonFACE(f.flag = 1, constraintDF_FACE, equilDF_FACE)

#### Step 2 store run 1 - 10 constrainDF dataframe

### Run 1
constraintDF_FACE <- Perform_Analytical_Run1_AmazonFACE(f.flag = 2, constraintDF_FACE, equilDF_FACE)

### Run 2
constraintDF_FACE <- Perform_Analytical_Run2_AmazonFACE(f.flag = 2, constraintDF_FACE, equilDF_FACE)



#### Step 3 store run 1 - 10 equilDF dataframes

### Run 1
equilDF_FACE <- Perform_Analytical_Run1_AmazonFACE(f.flag = 3, constraintDF_FACE, equilDF_FACE)

### Run 2
equilDF_FACE <- Perform_Analytical_Run2_AmazonFACE(f.flag = 3, constraintDF_FACE, equilDF_FACE)


