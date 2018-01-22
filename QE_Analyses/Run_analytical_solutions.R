
#### This is the main script to call individual analytical run scripts and functions
#### under the QE framework
####
#### Run definitions:
#### Run 1: baseline, variable wood NC, baseline N cycle model
#### Run 2: baseline, fixed wood NC, baseline N cycle model

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

#### Step 2: return data list
r1 <- Perform_Analytical_Run1(f.flag = 2)
r2 <- Perform_Analytical_Run2(f.flag = 2)

################################################################################
#### End


