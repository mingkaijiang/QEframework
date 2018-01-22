
#### This is the main script to call individual analytical run scripts and functions
#### under the QE framework
####
#### Run definitions:
#### Run 1: baseline, variable wood NC, baseline N cycle model
#### Run 2: baseline, fixed wood NC, baseline N cycle model
#### Run 3: baseline N cycle model, with explicit mineral N pool simulated, variable wood NC
#### Run 4: baseline N cycle model, with explicit mineral N pool simulated, fixed wood NC
#### Run 5: GDAY plant mineral N uptake function, variable wood NC
#### Run 6: GDAY plant mineral N uptake function, fixed wood NC
#### Run 7: O-CN plant mineral N uptake function, variable wood NC
#### Run 8: O-CN plant mineral N uptake function, fixed wood NC

################################################################################
#### Source QE scripts

# There is a line in prepare_R that repeats this. For now, leave it here.
scriptfiles <- dir("QE_Scripts", pattern="[.]R$", recursive = TRUE, full.names = TRUE)
for(b in scriptfiles)source(b)

################################################################################
### f.flag: = 1 return saved plots
###         = 2 return list of two dataframes

#### Step 1: return plots
Perform_Analytical_Run1(f.flag = 1)
Perform_Analytical_Run2(f.flag = 1)
Perform_Analytical_Run3(f.flag = 1)
Perform_Analytical_Run4(f.flag = 1)
Perform_Analytical_Run5(f.flag = 1)
Perform_Analytical_Run6(f.flag = 1)
Perform_Analytical_Run7(f.flag = 1)
Perform_Analytical_Run8(f.flag = 1)

#### Step 2: return data list
r1 <- Perform_Analytical_Run1(f.flag = 2)
r2 <- Perform_Analytical_Run2(f.flag = 2)
r3 <- Perform_Analytical_Run3(f.flag = 2)
r4 <- Perform_Analytical_Run4(f.flag = 2)
r5 <- Perform_Analytical_Run5(f.flag = 2)
r6 <- Perform_Analytical_Run6(f.flag = 2)
r7 <- Perform_Analytical_Run7(f.flag = 2)
r8 <- Perform_Analytical_Run8(f.flag = 2)

################################################################################
#### End


