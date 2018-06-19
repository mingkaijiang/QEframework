
#### This is the main script to call individual analytical run scripts and functions
#### under the QE framework
####
#### Run definitions:
#### Run 1: baseline, variable wood NC, baseline N cycle model
#### Run 2: baseline, fixed wood NC, baseline N cycle model
#### Run 3: baseline N cycle model, with explicit mineral N pool simulated, variable wood NC
#### Run 4: baseline N cycle model, with explicit mineral N pool simulated, fixed wood NC
#### Run 5: GDAY plant mineral N uptake function, variable wood NC
#### Run 6: O-CN plant mineral N uptake function, variable wood NC
#### Run 7: Priming effect (i.e. exudation and faster turnover of slow SOM pool, fixed wood)
#### Run 8: Priming effect (i.e. exudation and faster turnover of slow SOM pool, variable wood)
#### Run 9: CLM baseline potential NPP, fixed wood NC
#### Run 10: CLM baseline potential NPP, variable wood NC
#### Run 11: CLM C cost of N uptake approach, fixed wood NC
#### Run 12: CLM C cost of N uptake approach, variable wood NC

################################################################################
### f.flag: = 1 return saved plots
###         = 2 return list of two dataframes

#### Step 1: return plots
#Perform_Analytical_Run1(f.flag = 1)          
#Perform_Analytical_Run2(f.flag = 1)
#Perform_Analytical_Run3(f.flag = 1)
#Perform_Analytical_Run31(f.flag = 1)
#Perform_Analytical_Run32(f.flag = 1)
#Perform_Analytical_Run33(f.flag = 1)
#Perform_Analytical_Run34(f.flag = 1)
#Perform_Analytical_Run4(f.flag = 1)
#Perform_Analytical_Run5(f.flag = 1)
#Perform_Analytical_Run51(f.flag = 1)  # varying ar
#Perform_Analytical_Run52(f.flag = 1)
#Perform_Analytical_Run53(f.flag = 1)  # varying sr
#Perform_Analytical_Run54(f.flag = 1)
#Perform_Analytical_Run6(f.flag = 1)
#Perform_Analytical_Run7(f.flag = 1)
#Perform_Analytical_Run70(f.flag = 1)
#Perform_Analytical_Run8(f.flag = 1)
#Perform_Analytical_Run80(f.flag = 1)
#Perform_Analytical_Run81(f.flag = 1)
#Perform_Analytical_Run82(f.flag = 1)
#Perform_Analytical_Run83(f.flag = 1)
#Perform_Analytical_Run84(f.flag = 1)
#Perform_Analytical_Run85(f.flag = 1)
#Perform_Analytical_Run86(f.flag = 1)
#Perform_Analytical_Run91(f.flag = 1)
#Perform_Analytical_Run101(f.flag = 1)
#Perform_Analytical_Run11(f.flag = 1)
#Perform_Analytical_Run12(f.flag = 1)
#Perform_Analytical_Run13(f.flag = 1)
#Perform_Analytical_Run15(f.flag = 1)
#Perform_Analytical_Run17(f.flag = 1)
#Perform_Analytical_Run19(f.flag = 1)

#### Step 2: return data list
r1 <- Perform_Analytical_Run1(f.flag = 2)
r2 <- Perform_Analytical_Run2(f.flag = 2)
r3 <- Perform_Analytical_Run3(f.flag = 2)
r31 <- Perform_Analytical_Run31(f.flag = 2)  # varying N uptake
r32 <- Perform_Analytical_Run32(f.flag = 2)
r33 <- Perform_Analytical_Run33(f.flag = 2)
r34 <- Perform_Analytical_Run34(f.flag = 2)
r4 <- Perform_Analytical_Run4(f.flag = 2)
r5 <- Perform_Analytical_Run5(f.flag = 2)
r51 <- Perform_Analytical_Run51(f.flag = 2)  # varying ar
r52 <- Perform_Analytical_Run52(f.flag = 2)
r53 <- Perform_Analytical_Run53(f.flag = 2)  # varying sr
r54 <- Perform_Analytical_Run54(f.flag = 2)
r6 <- Perform_Analytical_Run6(f.flag = 2)
r7 <- Perform_Analytical_Run7(f.flag = 2)
r70 <- Perform_Analytical_Run70(f.flag = 2)
r8 <- Perform_Analytical_Run8(f.flag = 2)
r80 <- Perform_Analytical_Run80(f.flag = 2)
r81 <- Perform_Analytical_Run81(f.flag = 2)
r82 <- Perform_Analytical_Run82(f.flag = 2)
r83 <- Perform_Analytical_Run83(f.flag = 2)
r84 <- Perform_Analytical_Run84(f.flag = 2)
r85 <- Perform_Analytical_Run85(f.flag = 2)
r86 <- Perform_Analytical_Run86(f.flag = 2)
r9 <- Perform_Analytical_Run91(f.flag = 2)
r10 <- Perform_Analytical_Run101(f.flag = 2)
#r11 <- Perform_Analytical_Run11(f.flag = 2)
#r12 <- Perform_Analytical_Run12(f.flag = 2)

################################################################################
#### End


