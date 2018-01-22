
#### This is the main script to call individual analytical run scripts and functions
#### under the QE framework
####
#### Run definitions:
#### Run 1: baseline, variable wood NC, baseline N cycle model
#### Run 2: baseline, fixed wood NC, baseline N cycle model

################################################################################
#### Create dataframes to store all the data
### for now, ignore this step as I haven't come up with a more efficient way of storing the outputs
### so, need to work on computing some outputs first, and see what storage strategy I should use thereafter. 

# constraintDF <- create_constraint_DF()
# equilDF <- create_equil_DF()

################################################################################
#### Source QE scripts

# There is a line in prepare_R that repeats this. For now, leave it here.
scriptfiles <- dir("QE_Scripts", pattern="[.]R$", recursive = TRUE, full.names = TRUE)
for(b in scriptfiles)source(b)

################################################################################
#### Step 1: simply run analytical solution and plot quasi-equil plots
### f.flag: = 1 simply plot analytical solution graph
###         = 2 return constraintDF
###         = 3 return equilDF
Perform_Analytical_Run1(f.flag = 1, constraintDF, equilDF)

################################################################################
#### Step 2 store run 1 - 10 constrainDF dataframe
### Run 1
constraintDF <- Perform_Analytical_Run1(f.flag = 2, constraintDF, equilDF)


################################################################################
#### Step 3 store run 1 - 10 equilDF dataframes
### Run 1
equilDF <- Perform_Analytical_Run1(f.flag = 3, constraintDF, equilDF)


################################################################################
#### End


