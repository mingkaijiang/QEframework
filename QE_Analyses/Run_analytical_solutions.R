
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
#### Run 3: same as Run1, but fixed wood stoichiometry
#### Run 4: same as Run1, but autotrophic respiration as a function of plant N concentration
#### Run 5: same as Run1, but with simplified LUE model
#### Run 6: same as Run2, but with simplfieid LUE model
#### Run 7: same as Run1, but with explicit mineral N pool
#### Run 8.1: same as Run7, but with N only model, nuptake ~ root biomass - GDAY approach
#### Run 8.2, same as Run8.1, but with N only model, O-CN approach (i.e. saturating function of mineral N)
#### Run 9: same as Run 7, but with N only model, fixed wood stoichiometry, and priming on
################################################################################
#### Create dataframes to store all the data
constraintDF <- create_constraint_DF()
equilDF <- create_equil_DF()

#### Step 1: simply run analytical solution and plot quasi-equil plots
### f.flag: = 1 simply plot analytical solution graph
###         = 2 return constraintDF
###         = 3 return equilDF
Perform_Analytical_Run1(f.flag = 1, constraintDF, equilDF)
Perform_Analytical_Run2(f.flag = 1, constraintDF, equilDF)
Perform_Analytical_Run3(f.flag = 1, constraintDF, equilDF)
Perform_Analytical_Run4(f.flag = 1, constraintDF, equilDF)
Perform_Analytical_Run5(f.flag = 1, constraintDF, equilDF)
Perform_Analytical_Run6(f.flag = 1, constraintDF, equilDF)
Perform_Analytical_Run7(f.flag = 1, constraintDF, equilDF)
Perform_Analytical_Run8_1(f.flag = 1, constraintDF, equilDF)
Perform_Analytical_Run8_2()
Perform_Analytical_Run9(f.flag = 1, constraintDF, equilDF)
Perform_Analytical_Run10(f.flag = 1, constraintDF, equilDF)
Perform_Analytical_Run11(f.flag = 1, constraintDF, equilDF)

#### Step 2 store run 1 - 10 constrainDF dataframe

### Run 1
constraintDF <- Perform_Analytical_Run1(f.flag = 2, constraintDF, equilDF)

### Run 2
constraintDF <- Perform_Analytical_Run2(f.flag = 2, constraintDF, equilDF)

### Run 3
constraintDF <- Perform_Analytical_Run3(f.flag = 2, constraintDF, equilDF)

#### Run 4
#constraintDF <- Perform_Analytical_Run4(f.flag = 2, constraintDF, equilDF)
#
#### Run 5
#constraintDF <- Perform_Analytical_Run5(f.flag = 2, constraintDF, equilDF)
#
#### Run 6
#constraintDF <- Perform_Analytical_Run6(f.flag = 2, constraintDF, equilDF)
#
#### Run 7
#constraintDF <- Perform_Analytical_Run7(f.flag = 2, constraintDF, equilDF)
#
#### Run 8
#constraintDF <- Perform_Analytical_Run8_1(f.flag = 2, constraintDF, equilDF)
#
#### Run 9
#constraintDF <- Perform_Analytical_Run9(f.flag = 2, constraintDF, equilDF)

#### Step 3 store run 1 - 10 equilDF dataframes

### Run 1
equilDF <- Perform_Analytical_Run1(f.flag = 3, constraintDF, equilDF)

### Run 2
equilDF <- Perform_Analytical_Run2(f.flag = 3, constraintDF, equilDF)

### Run 3
equilDF <- Perform_Analytical_Run3(f.flag = 3, constraintDF, equilDF)

#### Run 4
#equilDF <- Perform_Analytical_Run4(f.flag = 3, constraintDF, equilDF)
#
#### Run 5
#equilDF <- Perform_Analytical_Run5(f.flag = 3, constraintDF, equilDF)
#
#### Run 6
#equilDF <- Perform_Analytical_Run6(f.flag = 3, constraintDF, equilDF)
#
#### Run 7
#equilDF <- Perform_Analytical_Run7(f.flag = 3, constraintDF, equilDF)
#
#### Run 8
#equilDF <- Perform_Analytical_Run8_1(f.flag = 3, constraintDF, equilDF)
#
#### Run 9
#equilDF <- Perform_Analytical_Run9(f.flag = 3, constraintDF, equilDF)


