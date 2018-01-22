
#### Functions to generate Table 1
#### Purpose:
#### To generate % diff between aCO2 and eCO2 over different timescales
#### only comparing N and NP models (i.e. Analytical run 1 and 2)
####
################################################################################

######################## Functions ###################################
compute_Table_1 <- function(destDir) {
    
    ######## Make the output table
    out.tab <- matrix(ncol=7, nrow = 8)
    colnames(out.tab) <- c("model", "NPP_350", "NPP_700", "I", "M", "L", "VL")
    out.tab <- as.data.frame(out.tab)
    out.tab[1,"model"] <- "Baseline N-P model"
    out.tab[2,"model"] <- "N-limited model"
    out.tab[3,"model"] <- "Fixed wood"
    out.tab[4,"model"] <- "Explicit N uptake, fixed coefficient"
    out.tab[5,"model"] <- "Explicit N uptake, GDAY"
    out.tab[6,"model"] <- "Explicit N uptake, OCN"
    out.tab[7,"model"] <- "Exudation on, priming off"
    out.tab[8,"model"] <- "Exudation on, priming on"
    
    
    ######### Run analytical run 1
    source("Parameters/Analytical_Run1_Parameters.R")
    
    # create a range of nc for shoot to initiate
    nfseq <- round(seq(0.001, 0.1, by = 0.001),5)
    a_nf <- as.data.frame(allocn(nfseq))
    
    # using very long term relationship to calculate pf from nf
    pfseq <- inferpfVL(nfseq, a_nf)
    a_pf <- as.data.frame(allocp(pfseq))
    
    # calculate photosynthetic constraint at CO2 = 350
    Photo350 <- photo_constraint_full_cnp(nfseq, pfseq, a_nf, a_pf, CO2_1)
    
    ### calculate very long term NC and PC constraint on NPP, respectively
    NCVLONG <- VLong_constraint_N(nf=nfseq, nfdf=a_nf)
    
    ### NPP derived from PCVLONG should match NPP from NCVLONG
    PCVLONG <- VLong_constraint_P(pf=pfseq, pfdf=a_pf)
    
    ### finding the equilibrium point between photosynthesis and very long term nutrient constraints
    VLong_equil <- solveVLong_full_cnp(CO2=CO2_1)
    
    ### Get Cpassive from very-long nutrient cycling solution
    aequiln <- allocn(VLong_equil$equilnf)
    aequilp <- allocp(VLong_equil$equilpf)
    #pass <- passive(df=VLong_equil$equilnf, a=aequiln)
    #omega <- aequiln$af*pass$omegaf + aequiln$ar*pass$omegar
    #CpassVLong <- omega*VLong_equil$equilNPP/pass$decomp/(1-pass$qq)*1000.0
    
    pass <- slow_pool(df=VLong_equil$equilnf, a=aequiln)
    omegap <- aequiln$af*pass$omegafp + aequiln$ar*pass$omegarp
    CpassVLong <- omegap*VLong_equil$equilNPP/pass$decomp_p/(1-pass$qpq)*1000.0
    
    # Calculate long term nutrient constraint
    NCLONG <- Long_constraint_N(nfseq, a_nf, CpassVLong,
                                NinL = Nin)#+NrelwoodVLong)
    
    # Calculate pf based on nf of long-term nutrient exchange
    pfseqL <- inferpfL(nfseq, a_nf, PinL = Pin,#+PrelwoodVLong,
                       NinL = Nin,#+NrelwoodVLong,
                       Cpass=CpassVLong)
    
    PCLONG <- Long_constraint_P(nfseq, pfseqL, allocp(pfseqL),
                                CpassVLong, PinL=Pin)#+PrelwoodVLong)
    
    # Find long term equilibrium point
    Long_equil <- solveLong_full_cnp(CO2=CO2_1, Cpass=CpassVLong, NinL = Nin,#+NrelwoodVLong, 
                                     PinL=Pin)#+PrelwoodVLong)
    
    # Get Cslow from long nutrient cycling solution
    omegas <- aequiln$af*pass$omegafs + aequiln$ar*pass$omegars
    CslowLong <- omegas*Long_equil$equilNPP/pass$decomp_s/(1-pass$qsq)*1000.0
    
    ### Calculate nutrient release from slow woody pool
    PrelwoodVLong <- aequilp$aw*aequilp$pw*VLong_equil$equilNPP*1000.0
    NrelwoodVLong <- aequiln$aw*aequiln$nw*VLong_equil$equilNPP*1000.0
    
    # Calculate pf based on nf of medium-term nutrient exchange
    pfseqM <- inferpfM(nfseq, a_nf, PinM = Pin+PrelwoodVLong,
                       NinM = Nin+NrelwoodVLong,
                       CpassL=CpassVLong, CpassM=CslowLong)
    
    # Calculate medium term nutrient constraint
    NCMEDIUM <- NConsMedium(df=nfseq, 
                            a=a_nf, 
                            Cpass=CpassVLong, 
                            Cslow=CslowLong, 
                            NinL = Nin+NrelwoodVLong)
    
    Medium_equil_350 <- solveMedium_full_cnp(CO2=CO2_1, Cpass=CpassVLong, Cslow=CslowLong, NinL = Nin+NrelwoodVLong,
                                             PinL=Pin+PrelwoodVLong)

    out350DF <- data.frame(nfseq, pfseq, pfseqL, Photo350, NCVLONG, NCLONG)
    colnames(out350DF) <- c("nc", "pc_VL", "pc_350_L", "NPP_350", "NPP_VL",
                            "nleach_VL", "NPP_350_L", "nwood_L", "nburial_L",
                            "nleach_L", "aw")
    equil350DF <- data.frame(VLong_equil, Long_equil)
    colnames(equil350DF) <- c("nc_VL", "pc_VL", "NPP_VL", 
                              "nc_L","pc_L", "NPP_L")
    
    ##### CO2 = 700
    
    # N:C and P:C ratio
    nfseq <- round(seq(0.001, 0.1, by = 0.001),5)
    a_nf <- as.data.frame(allocn(nfseq))
    
    # using very long term relationship to calculate pf from nf
    pfseq <- inferpfVL(nfseq, a_nf)
    a_pf <- as.data.frame(allocp(pfseq))
    
    # calculate NC vs. NPP at CO2 = 350 respectively
    Photo700 <- photo_constraint_full_cnp(nfseq, pfseq, a_nf, a_pf, CO2_2)
    
    ### calculate very long term NC and PC constraint on NPP, respectively
    NCVLONG <- VLong_constraint_N(nf=nfseq, nfdf=a_nf)
    
    ### NPP derived from PCVLONG should match NPP from NCVLONG
    PCVLONG <- VLong_constraint_P(pf=pfseq, pfdf=a_pf)
    
    ### finding the equilibrium point between photosynthesis and very long term nutrient constraints
    VLong_equil <- solveVLong_full_cnp(CO2=CO2_2)
    
    # Find long term equilibrium point
    Long_equil <- solveLong_full_cnp(CO2=CO2_2, Cpass=CpassVLong, NinL = Nin,#+NrelwoodVLong, 
                                     PinL=Pin)#+PrelwoodVLong)
    
    # Find medium term equilibrium point
    Medium_equil_350 <- solveMedium_full_cnp(CO2_1, Cpass = CpassVLong, Cslow = CslowLong, 
                                             NinL=Nin+NrelwoodVLong, PinL=Pin+PrelwoodVLong)
    Medium_equil_700 <- solveMedium_full_cnp(CO2_2, Cpass = CpassVLong, Cslow = CslowLong, 
                                             NinL=Nin+NrelwoodVLong, PinL=Pin+PrelwoodVLong)
    
    out700DF <- data.frame(nfseq, pfseq, pfseqL, Photo700, NCVLONG, NCLONG)
    colnames(out700DF) <- c("nc", "pc_VL", "pc_700_L", "NPP_700", "NPP_VL",
                            "nleach_VL", "NPP_700_L", "nwood_L", "nburial_L",
                            "nleach_L", "aw")
    
    equil700DF <- data.frame(VLong_equil, Long_equil)
    colnames(equil700DF) <- c("nc_VL", "pc_VL", "NPP_VL", 
                              "nc_L","pc_L", "NPP_L")
    
    # get the point instantaneous NPP response to doubling of CO2
    df700 <- as.data.frame(cbind(round(nfseq,3), Photo700))
    inst700 <- inst_NPP(equil350DF$nc_VL, df700)
    
    ### pass on information to the output data frame
    out.tab[1,"NPP_350"] <- equil350DF$NPP_VL
    out.tab[1,"NPP_700"] <- equil700DF$NPP_VL
    out.tab[1,"I"] <- (inst700$equilNPP - equil350DF$NPP_VL) / equil350DF$NPP_VL
    out.tab[1,"M"] <- (Medium_equil_700$equilNPP - equil350DF$NPP_VL) / equil350DF$NPP_VL
    out.tab[1,"L"] <- (equil700DF$NPP_L - equil350DF$NPP_VL) / equil350DF$NPP_VL
    out.tab[1,"VL"] <- (equil700DF$NPP_VL - equil350DF$NPP_VL) / equil350DF$NPP_VL
    
    ######### Run analytical run 2
    source("Parameters/Analytical_Run2_Parameters.R")
    
    # N:C ratios for x-axis
    nfseq <- seq(0.001,0.1,by=0.001)
    # need allocation fractions here
    a_vec <- allocn(nfseq)
    
    # plot photosynthetic constraints
    PC350 <- photo_constraint_full_cn(nfseq,a_vec,CO2=CO2_1)
    PC700 <- photo_constraint_full_cn(nfseq,a_vec,CO2=CO2_2)
    
    #plot very-long nutrient cycling constraint
    NCVLONG <- VLong_constraint_N(nfseq,a_vec)
    
    #solve very-long nutrient cycling constraint
    VLong <- solveVLong_full_cn(CO2=CO2_1)
    #get Cpassive from very-long nutrient cycling solution
    aequil <- allocn(VLong$equilnf)
    pass <- slow_pool(df=VLong$equilnf, a=aequil)
    omegap <- aequil$af*pass$omegafp + aequil$ar*pass$omegarp
    omegas <- aequil$af*pass$omegafs + aequil$ar*pass$omegars
    
    CpassVLong <- omegap*VLong$equilNPP/pass$decomp_p/(1-pass$qpq)*1000.0
    
    NrelwoodVLong <- aequil$aw*aequil$nw*VLong$equilNPP*1000
    
    #now plot long-term constraint with this Cpassive
    NCHUGH <- Long_constraint_N(nfseq,a_vec,Cpass=CpassVLong,Nin)#+NrelwoodVLong)
    
    # Solve longterm equilibrium
    equil_long_350 <- solveLong_full_cn(CO2=CO2_1, Cpass=CpassVLong, NinL = Nin)#+NrelwoodVLong)
    equil_long_700 <- solveLong_full_cn(CO2=CO2_2, Cpass=CpassVLong, NinL = Nin)#+NrelwoodVLong)
    
    CslowLong <- omegas*equil_long_350$equilNPP/pass$decomp_s/(1-pass$qsq)*1000.0
    
    # plot medium nutrient cycling constraint
    NCMEDIUM <- NConsMedium(nfseq, a_vec, CpassVLong, CslowLong, Nin+NrelwoodVLong)
    
    # solve medium term equilibrium at CO2 = 700 ppm
    equil_medium_700 <- solveMedium(CO2_2,Cpass=CpassVLong,Cslow=CslowLong,Nin=Nin+NrelwoodVLong)
    
    # get the point instantaneous NPP response to doubling of CO2
    df700 <- as.data.frame(cbind(round(nfseq,3), PC700))
    inst700 <- inst_NPP(VLong$equilnf, df700)
    
    ## locate the intersect between VL nutrient constraint and CO2 = 700
    VLong700 <- solveVLong_full_cn(CO2=CO2_2)
    
    ### pass on information to the output data frame
    out.tab[2,"NPP_350"] <- VLong$equilNPP
    out.tab[2,"NPP_700"] <- VLong700$equilNPP
    out.tab[2,"I"] <- (inst700$equilNPP - VLong$equilNPP) / VLong$equilNPP
    out.tab[2,"M"] <- (equil_medium_700$equilNPP - VLong$equilNPP) / VLong$equilNPP
    out.tab[2,"L"] <- (equil_long_700$equilNPP - VLong$equilNPP) / VLong$equilNPP
    out.tab[2,"VL"] <- (VLong700$equilNPP - VLong$equilNPP) / VLong$equilNPP
    
    
    ######### Run analytical run 3
    source("Parameters/Analytical_Run3_Parameters.R")
    
    # create a range of nc for shoot to initiate
    nfseq <- round(seq(0.001, 0.1, by = 0.001),5)
    a_nf <- as.data.frame(allocn(nfseq))
    
    # using very long term relationship to calculate pf from nf
    pfseq <- inferpfVL(nfseq, a_nf)
    a_pf <- as.data.frame(allocp(pfseq))
    
    # calculate photosynthetic constraint at CO2 = 350
    Photo350 <- photo_constraint_full_cnp(nfseq, pfseq, a_nf, a_pf, CO2_1)
    
    ### calculate very long term NC and PC constraint on NPP, respectively
    NCVLONG <- VLong_constraint_N(nf=nfseq, nfdf=a_nf)
    
    ### NPP derived from PCVLONG should match NPP from NCVLONG
    PCVLONG <- VLong_constraint_P(pf=pfseq, pfdf=a_pf)
    
    ### finding the equilibrium point between photosynthesis and very long term nutrient constraints
    VLong_equil <- solveVLong_full_cnp(CO2=CO2_1)
    
    ### Get Cpassive from very-long nutrient cycling solution
    aequiln <- allocn(VLong_equil$equilnf)
    aequilp <- allocp(VLong_equil$equilpf)
    
    pass <- slow_pool(df=VLong_equil$equilnf, a=aequiln)
    omegap <- aequiln$af*pass$omegafp + aequiln$ar*pass$omegarp
    CpassVLong <- omegap*VLong_equil$equilNPP/pass$decomp_p/(1-pass$qpq)*1000.0
    
    # Calculate long term nutrient constraint
    NCLONG <- Long_constraint_N(nfseq, a_nf, CpassVLong,
                                NinL = Nin)#+NrelwoodVLong)
    
    # Calculate pf based on nf of long-term nutrient exchange
    pfseqL <- inferpfL(nfseq, a_nf, PinL = Pin,#+PrelwoodVLong,
                       NinL = Nin,#+NrelwoodVLong,
                       Cpass=CpassVLong)
    
    PCLONG <- Long_constraint_P(nfseq, pfseqL, allocp(pfseqL),
                                CpassVLong, PinL=Pin)#+PrelwoodVLong)
    
    # Find long term equilibrium point
    Long_equil <- solveLong_full_cnp(CO2=CO2_1, Cpass=CpassVLong, NinL = Nin,#+NrelwoodVLong, 
                                     PinL=Pin)#+PrelwoodVLong)
    
    # Get Cslow from long nutrient cycling solution
    omegas <- aequiln$af*pass$omegafs + aequiln$ar*pass$omegars
    CslowLong <- omegas*Long_equil$equilNPP/pass$decomp_s/(1-pass$qsq)*1000.0
    
    ### Calculate nutrient release from recalcitrant pools
    PrelwoodVLong <- aequilp$aw*aequilp$pw*VLong_equil$equilNPP*1000.0
    NrelwoodVLong <- aequiln$aw*aequiln$nw*VLong_equil$equilNPP*1000.0
    
    # Calculate pf based on nf of medium-term nutrient exchange
    pfseqM <- inferpfM(nfseq, a_nf, PinM = Pin+PrelwoodVLong,
                       NinM = Nin+NrelwoodVLong,
                       CpassL=CpassVLong, CpassM=CslowLong)
    
    # Calculate medium term nutrient constraint
    NCMEDIUM <- NConsMedium(df=nfseq, 
                            a=a_nf, 
                            Cpass=CpassVLong, 
                            Cslow=CslowLong, 
                            NinL = Nin+NrelwoodVLong)
    
    Medium_equil_350 <- solveMedium_full_cnp(CO2=CO2_1, Cpass=CpassVLong, Cslow=CslowLong, NinL = Nin+NrelwoodVLong,
                                             PinL=Pin+PrelwoodVLong)
    
    out350DF <- data.frame(nfseq, pfseq, pfseqL, Photo350, NCVLONG, NCLONG)
    colnames(out350DF) <- c("nc", "pc_VL", "pc_350_L", "NPP_350", "NPP_VL",
                            "nleach_VL", "NPP_350_L", "nwood_L", "nburial_L",
                            "nleach_L", "aw")
    equil350DF <- data.frame(VLong_equil, Long_equil)
    colnames(equil350DF) <- c("nc_VL", "pc_VL", "NPP_VL", 
                              "nc_L","pc_L", "NPP_L")
    
    ##### CO2 = 700
    nfseq <- round(seq(0.001, 0.1, by = 0.001),5)
    a_nf <- as.data.frame(allocn(nfseq))
    
    # using very long term relationship to calculate pf from nf
    pfseq <- inferpfVL(nfseq, a_nf)
    a_pf <- as.data.frame(allocp(pfseq))
    
    # calculate NC vs. NPP at CO2 = 350 respectively
    Photo700 <- photo_constraint_full_cnp(nfseq, pfseq, a_nf, a_pf, CO2_2)
    
    ### calculate very long term NC and PC constraint on NPP, respectively
    NCVLONG <- VLong_constraint_N(nf=nfseq, nfdf=a_nf)
    
    ### NPP derived from PCVLONG should match NPP from NCVLONG
    PCVLONG <- VLong_constraint_P(pf=pfseq, pfdf=a_pf)
    
    ### finding the equilibrium point between photosynthesis and very long term nutrient constraints
    VLong700 <- solveVLong_full_cnp(CO2=CO2_2)
    
    # Find long term equilibrium point
    Long_equil <- solveLong_full_cnp(CO2=CO2_2, Cpass=CpassVLong, NinL = Nin,#+NrelwoodVLong, 
                                     PinL=Pin)#+PrelwoodVLong)
    
    # Find medium term equilibrium point
    Medium_equil_350 <- solveMedium_full_cnp(CO2_1, Cpass = CpassVLong, Cslow = CslowLong, 
                                             NinL=Nin+NrelwoodVLong, PinL=Pin+PrelwoodVLong)
    Medium_equil_700 <- solveMedium_full_cnp(CO2_2, Cpass = CpassVLong, Cslow = CslowLong, 
                                             NinL=Nin+NrelwoodVLong, PinL=Pin+PrelwoodVLong)
    
    out700DF <- data.frame(nfseq, pfseq, pfseqL, Photo700, NCVLONG, NCLONG)
    colnames(out700DF) <- c("nc", "pc_VL", "pc_700_L", "NPP_700", "NPP_VL",
                            "nleach_VL", "NPP_700_L", "nwood_L", "nburial_L",
                            "nleach_L", "aw")
    
    equil700DF <- data.frame(VLong700, Long_equil)
    colnames(equil700DF) <- c("nc_VL", "pc_VL", "NPP_VL", 
                              "nc_L","pc_L", "NPP_L")
    
    # get the point instantaneous NPP response to doubling of CO2
    df700 <- as.data.frame(cbind(round(nfseq,3), Photo700))
    inst700 <- inst_NPP(equil350DF$nc_VL, df700)
    
    out.tab[3,"NPP_350"] <- equil350DF$NPP_VL
    out.tab[3,"NPP_700"] <- equil700DF$NPP_VL
    out.tab[3,"I"] <- (inst700$equilNPP - equil350DF$NPP_VL) / equil350DF$NPP_VL
    out.tab[3,"M"] <- (Medium_equil_700$equilNPP - equil350DF$NPP_VL) / equil350DF$NPP_VL
    out.tab[3,"L"] <- (equil700DF$NPP_L - equil350DF$NPP_VL) / equil350DF$NPP_VL
    out.tab[3,"VL"] <- (equil700DF$NPP_VL - equil350DF$NPP_VL) / equil350DF$NPP_VL
    
    ######## Run analytical run 7 with varying nrateuptake coefficients - 1
    source("Parameters/Analytical_Run7_Parameters.R")
    
    # create nc and pc for shoot to initiate
    nfseq <- round(seq(0.001, 0.1, by = 0.001),5)
    a_nf <- as.data.frame(allocn(nfseq))
    
    #### Uptake rate = baseline = 1
    
    ##### CO2 = 350
    # calculate NC vs. NPP at CO2 = 350 respectively
    NC350 <- photo_constraint_full_cn(nfseq, a_nf, CO2_1)
    
    # calculate very long term NC and PC constraint on NPP, respectively
    NCVLONG <- NConsVLong_expl_min(df=nfseq,a=a_nf)
    
    # solve very-long nutrient cycling constraint
    VLongN <- solveVLong_expl_min(CO2_1)
    
    # Get Cpassive from very-long nutrient cycling solution
    aequiln <- allocn(VLongN$equilnf)
    
    pass <- slow_pool(df=VLongN$equilnf, a=aequiln)
    omegap <- aequiln$af*pass$omegafp + aequiln$ar*pass$omegarp
    CpassVLong <- omegap*VLongN$equilNPP/pass$decomp_p/(1-pass$qpq)*1000.0
    
    # Calculate long term nutrieng constraint
    NCHUGH <- NConsLong_expl_min(nfseq, a_nf,CpassVLong,
                                 NinL = Nin)#+NrelwoodVLong)
    
    # Find equilibrate intersection and plot
    LongN <- solveLong_expl_min(CO2_1, Cpass=CpassVLong, NinL= Nin)#+PrelwoodVLong)
    
    # Get Cslow from long nutrient cycling solution
    omegas <- aequiln$af*pass$omegafs + aequiln$ar*pass$omegars
    CslowLong <- omegas*LongN$equilNPP/pass$decomp_s/(1-pass$qsq)*1000.0
    
    # Calculate nutrient release from recalcitrant pools
    NrelwoodVLong <- aequiln$aw*aequiln$nw*VLongN$equilNPP_N*1000.0
    
    # Calculate medium term nutrieng constraint
    NCMEDIUM_1 <- NConsMedium_expl_min(nfseq, a_nf,CpassVLong, CslowLong,
                                       NinL = Nin+NrelwoodVLong)
    
    out350DF_1 <- data.frame(nfseq, NC350, NCVLONG, NCHUGH)
    colnames(out350DF_1) <- c("nc", "NPP_350", "NPP_VL",
                              "nleach_VL", "NPP_350_L", "nwood_L", "nburial_L",
                              "nleach_L", "aw")
    equil350DF_1 <- data.frame(VLongN, LongN)
    colnames(equil350DF_1) <- c("nc_VL", "NPP_VL", 
                                "nc_L", "NPP_L")
    
    ##### CO2 = 700
    # calculate NC vs. NPP at CO2 = 350 respectively
    NC700 <- photo_constraint_full_cn(nfseq, a_nf, CO2_2)
    
    # calculate very long term NC and PC constraint on NPP, respectively
    NCVLONG <- NConsVLong_expl_min(df=nfseq,a=a_nf)
    
    # solve very-long nutrient cycling constraint
    VLongN <- solveVLong_expl_min(CO2_2)
    
    out700DF_1 <- data.frame(nfseq, NC700, NCVLONG, NCHUGH)
    colnames(out700DF_1) <- c("nc", "NPP_700", "NPP_VL",
                              "nleach_VL", "NPP_700_L", "nwood_L", "nburial_L",
                              "nleach_L", "aw")
    
    # Find equilibrate intersection and plot
    LongN <- solveLong_expl_min(CO2_2, Cpass=CpassVLong, NinL=Nin)#+NrelwoodVLong)
    
    equil700DF_1 <- data.frame(VLongN, LongN)
    colnames(equil700DF_1) <- c("nc_VL", "NPP_VL", 
                                "nc_L", "NPP_L")
    
    # Find medium term equilibrium point
    Medium_equil_350_1 <- solveMedium_expl_min(CO2_1, Cpass = CpassVLong, Cslow = CslowLong, 
                                               NinL=Nin+NrelwoodVLong)
    Medium_equil_700_1 <- solveMedium_expl_min(CO2_2, Cpass = CpassVLong, Cslow = CslowLong, 
                                               NinL=Nin+NrelwoodVLong)
    
    # get the point instantaneous NPP response to doubling of CO2
    df700 <- as.data.frame(cbind(round(nfseq,3), NC700))
    inst700_1 <- inst_NPP(equil350DF_1$nc_VL, df700)
    
    out.tab[4,"NPP_350"] <- equil350DF_1$NPP_VL
    out.tab[4,"NPP_700"] <- equil700DF_1$NPP_VL
    out.tab[4,"I"] <- (inst700_1$equilNPP - equil350DF_1$NPP_VL) / equil350DF_1$NPP_VL
    out.tab[4,"M"] <- (Medium_equil_700_1$equilNPP - equil350DF_1$NPP_VL) / equil350DF_1$NPP_VL
    out.tab[4,"L"] <- (equil700DF_1$NPP_L - equil350DF_1$NPP_VL) / equil350DF_1$NPP_VL
    out.tab[4,"VL"] <- (equil700DF_1$NPP_VL - equil350DF_1$NPP_VL) / equil350DF_1$NPP_VL
    

    ######## Run analytical run 8.1 and 8.2 with GDAY and OCN approaches
    ######## GDAY approach
    source("Parameters/Analytical_Run8_1_Parameters.R")
    
    # N:C ratios for x-axis
    nfseq <- seq(0.001,0.1,by=0.001)
    # need allocation fractions here
    a_vec <- allocn(nfseq)
    
    # plot photosynthetic constraints
    PC350_gday <- photo_constraint_full_cn(nfseq,a_vec,CO2=CO2_1)
    PC700_gday <- photo_constraint_full_cn(nfseq,a_vec,CO2=CO2_2)
    
    # calculate very long term NC and PC constraint on NPP, respectively
    NCVLONG_gday <- NConsVLong_root_gday(df=nfseq,a=a_vec)
    
    # solve very-long nutrient cycling constraint
    VLongN_gday <- solveVLong_root_gday(CO2_1)
    
    # Get Cpassive from very-long nutrient cycling solution
    aequiln <- allocn(VLongN_gday$equilnf)
    pass <- slow_pool(df=VLongN_gday$equilnf, a=aequiln)
    omegap <- aequiln$af*pass$omegafp + aequiln$ar*pass$omegarp
    omegas <- aequiln$af*pass$omegafs + aequiln$ar*pass$omegars
    CpassVLong <- omegap*VLongN_gday$equilNPP/pass$decomp_p/(1-pass$qpq)*1000.0
    
    # Calculate nutrient release from recalcitrant pools
    NrelwoodVLong <- aequiln$aw*aequiln$nw*VLongN_gday$equilNPP*1000.0
    
    # Calculate long term nutrieng constraint
    NCHUGH_gday <- NConsLong_root_gday(df=nfseq, a=a_vec,Cpass=CpassVLong,
                                       NinL = Nin)#+NrelwoodVLong)
    
    # Find equilibrate intersection and plot
    equil_long_350_gday <- solveLong_root_gday(CO2_1, Cpass=CpassVLong, NinL= Nin)#+NrelwoodVLong)
    equil_long_700_gday <- solveLong_root_gday(CO2_2, Cpass=CpassVLong, NinL= Nin)#+NrelwoodVLong)
    
    CslowLong <- omegas*equil_long_350_gday$equilNPP/pass$decomp_s/(1-pass$qsq)*1000.0
    
    # plot medium nutrient cycling constraint
    NCMEDIUM_gday <- NConsMedium_root_gday(nfseq, a_vec, Cpass=CpassVLong, Cslow=CslowLong, NinL=Nin+NrelwoodVLong)
    
    # solve medium term equilibrium at CO2 = 700 ppm
    equil_medium_700_gday <- solveMedium_root_gday(CO2_2,Cpass=CpassVLong,Cslow=CslowLong,Nin=Nin+NrelwoodVLong)
    
    # get the point instantaneous NPP response to doubling of CO2
    df700 <- as.data.frame(cbind(round(nfseq,3), PC700_gday))
    inst700_gday <- inst_NPP(VLongN_gday$equilnf, df700)
    
    ## locate the intersect between VL nutrient constraint and CO2 = 700
    VLong700_gday <- solveVLong_root_gday(CO2=CO2_2)
    
    out.tab[5,"NPP_350"] <- VLongN_gday$equilNPP
    out.tab[5,"NPP_700"] <- VLong700_gday$equilNPP
    out.tab[5,"I"] <- (inst700_gday$equilNPP - VLongN_gday$equilNPP) / VLongN_gday$equilNPP
    out.tab[5,"M"] <- (equil_medium_700_gday$equilNPP - VLongN_gday$equilNPP) / VLongN_gday$equilNPP
    out.tab[5,"L"] <- (equil_long_700_gday$equilNPP - VLongN_gday$equilNPP) / VLongN_gday$equilNPP
    out.tab[5,"VL"] <- ( VLong700_gday$equilNPP - VLongN_gday$equilNPP) / VLongN_gday$equilNPP
    
    ######## OCN approach
    source("Parameters/Analytical_Run8_2_Parameters.R")
    
    # N:C ratios for x-axis
    nfseq <- seq(0.001,0.1,by=0.001)
    # need allocation fractions here
    a_vec <- allocn(nfseq)
    
    # plot photosynthetic constraints
    PC350_ocn <- photo_constraint_full_cn(nfseq,a_vec,CO2=CO2_1)
    PC700_ocn <- photo_constraint_full_cn(nfseq,a_vec,CO2=CO2_2)
    
    # calculate very long term NC and PC constraint on NPP, respectively
    VLongN_ocn <- NConsVLong_root_ocn(CO2_1)
    
    # Get Cpassive from very-long nutrient cycling solution
    aequiln <- allocn(VLongN_ocn$equilnf)
    pass <- slow_pool(df=VLongN_ocn$equilnf, a=aequiln)
    omegap <- aequiln$af*pass$omegafp + aequiln$ar*pass$omegarp
    omegas <- aequiln$af*pass$omegafs + aequiln$ar*pass$omegars
    CpassVLong <- omegap*VLongN_ocn$equilNPP/pass$decomp_p/(1-pass$qpq)*1000.0
    
    # Calculate nutrient release from recalcitrant pools
    NrelwoodVLong <- aequiln$aw*aequiln$nw*VLongN_ocn$equilNPP*1000.0
    
    # Calculate long term nutrieng constraint
    NCHUGH_ocn <- NConsLong_root_ocn(df=nfseq, a=a_vec,Cpass=CpassVLong,
                                     NinL = Nin)#+NrelwoodVLong)
    
    # Find equilibrate intersection and plot
    equil_long_350_ocn <- solveLong_root_ocn(CO2_1, Cpass=CpassVLong, NinL= Nin)#+NrelwoodVLong)
    equil_long_700_ocn <- solveLong_root_ocn(CO2_2, Cpass=CpassVLong, NinL= Nin)#+NrelwoodVLong)
    
    CslowLong <- omegas*equil_long_350_ocn$equilNPP/pass$decomp_s/(1-pass$qsq)*1000.0
    
    # plot medium nutrient cycling constraint
    NCMEDIUM_ocn <- NConsMedium_root_ocn(nfseq, a_vec, Cpass=CpassVLong, Cslow=CslowLong, NinL=Nin+NrelwoodVLong)
    
    # solve medium term equilibrium at CO2 = 700 ppm
    equil_medium_700_ocn <- solveMedium_root_ocn(CO2_2,Cpass=CpassVLong,Cslow=CslowLong,Nin=Nin+NrelwoodVLong)
    
    # get the point instantaneous NPP response to doubling of CO2
    df700 <- as.data.frame(cbind(round(nfseq,3), PC700_ocn))
    inst700_ocn <- inst_NPP(VLongN_ocn$equilnf, df700)
    
    ## locate the intersect between VL nutrient constraint and CO2 = 700
    VLong700_ocn <-  NConsVLong_root_ocn(CO2_2)
    
    out.tab[6,"NPP_350"] <- VLongN_ocn$equilNPP
    out.tab[6,"NPP_700"] <- VLong700_ocn$equilNPP
    out.tab[6,"I"] <- (inst700_ocn$equilNPP - VLongN_ocn$equilNPP) / VLongN_ocn$equilNPP
    out.tab[6,"M"] <- (equil_medium_700_ocn$equilNPP - VLongN_ocn$equilNPP) / VLongN_ocn$equilNPP
    out.tab[6,"L"] <- (equil_long_700_ocn$equilNPP - VLongN_ocn$equilNPP) / VLongN_ocn$equilNPP
    out.tab[6,"VL"] <- ( VLong700_ocn$equilNPP - VLongN_ocn$equilNPP) / VLongN_ocn$equilNPP
    
    
    ########  Analytical run 11 - exudation on, priming off
    source("Parameters/Analytical_Run10_Parameters.R")
    
    # create nc and pc for shoot to initiate
    nfseq_off <- round(seq(0.001, 0.1, by = 0.001),5)
    a_vec <- as.data.frame(allocn_exudation(nfseq_off))
    
    # plot photosynthetic constraints
    PC350_off <- photo_constraint_full_cn(nfseq_off,a_vec,CO2=CO2_1)
    PC700_off <- photo_constraint_full_cn(nfseq_off,a_vec,CO2=CO2_2)
    
    # calculate very long term NC and PC constraint on NPP, respectively
    NCVLONG_off <- NConsVLong(df=nfseq_off,a=a_vec)
    
    # solve very-long nutrient cycling constraint
    VLongN_off <- solveVLong_full_cn_medium(CO2_1)
    
    # Get Cpassive from very-long nutrient cycling solution
    aequiln <- allocn_exudation(VLongN_off$equilnf)
    pass <- slow_pool(df=VLongN_off$equilnf, a=aequiln)
    omegap <- aequiln$af*pass$omegafp + aequiln$ar*pass$omegarp
    CpassVLong <- omegap*VLongN_off$equilNPP/pass$decomp_p/(1-pass$qpq)*1000.0
    
    # Get Cslow from long nutrient cycling solution
    omegas <- aequiln$af*pass$omegafs + aequiln$ar*pass$omegars
    CslowLong <- omegas*VLongN_off$equilNPP/pass$decomp_s/(1-pass$qsq)*1000.0
    
    # Calculate nutrient release from woody pool
    NrelwoodVLong <- aequiln$aw*aequiln$nw*VLongN_off$equilNPP*1000.0
    
    # Calculate long term nutrient constraint
    NCHUGH_off <- NConsLong(nfseq_off, a_vec, CpassVLong,
                            NinL = Nin)
    
    # Calculate medium term nutrient constraint
    NCMEDIUM_off <- NConsMedium(df=nfseq_off, 
                                a=a_vec, 
                                Cpass=CpassVLong, 
                                Cslow=CslowLong, 
                                NinL = Nin+NrelwoodVLong)
    
    
    # Solve medium equilibrium
    equil_long_350_off <- solveLong_full_cn_medium(CO2=CO2_1, Cpass=CpassVLong, NinL = Nin)
    equil_long_700_off <- solveLong_full_cn_medium(CO2=CO2_2, Cpass=CpassVLong, NinL = Nin)
    
    # Solve medium equilibrium
    equil_medium_350_off <- solveMedium(CO2=CO2_1, Cpass=CpassVLong, Cslow=CslowLong, NinL = Nin+NrelwoodVLong)
    equil_medium_700_off <- solveMedium(CO2=CO2_2, Cpass=CpassVLong, Cslow=CslowLong, NinL = Nin+NrelwoodVLong)
    
    # get the point instantaneous NPP response to doubling of CO2
    df700 <- as.data.frame(cbind(round(nfseq_off,3), PC700_off))
    inst700_off <- inst_NPP(VLongN_off$equilnf, df700)
    
    ## locate the intersect between VL nutrient constraint and CO2 = 700
    VLong700_off <- solveVLong_full_cn_medium(CO2_2)
    
    out.tab[7,"NPP_350"] <- VLongN_off$equilNPP
    out.tab[7,"NPP_700"] <- VLong700_off$equilNPP
    out.tab[7,"I"] <- (inst700_off$equilNPP - VLongN_off$equilNPP) / VLongN_off$equilNPP
    out.tab[7,"M"] <- (equil_medium_700_off$equilNPP - VLongN_off$equilNPP) / VLongN_off$equilNPP
    out.tab[7,"L"] <- (equil_long_700_off$equilNPP - VLongN_off$equilNPP) / VLongN_off$equilNPP
    out.tab[7,"VL"] <- ( VLong700_off$equilNPP - VLongN_off$equilNPP) / VLongN_off$equilNPP

    
    
    ######## Analytical run 10 - exudation on, priming on
    # create nc and pc for shoot to initiate
    nfseq_on <- round(seq(0.001, 0.1, by = 0.001),5)
    a_vec <- as.data.frame(allocn_exudation(nfseq_on))
    
    # plot photosynthetic constraints
    PC350_on <- photo_constraint_full_cn(nfseq_on,a_vec,CO2=CO2_1)
    PC700_on <- photo_constraint_full_cn(nfseq_on,a_vec,CO2=CO2_2)
    
    # calculate very long term NC and PC constraint on NPP, respectively
    NCVLONG_on <- NConsVLong(df=nfseq_on,a=a_vec)
    
    # solve very-long nutrient cycling constraint
    VLongN_on <- solveVLong_full_cn_medium(CO2_1)
    
    # Get Cpassive from very-long nutrient cycling solution
    aequiln <- allocn_exudation(VLongN_on$equilnf)
    pass <- slow_pool(df=VLongN_on$equilnf, a=aequiln)
    omegap <- aequiln$af*pass$omegafp + aequiln$ar*pass$omegarp
    CpassVLong <- omegap*VLongN_on$equilNPP/pass$decomp_p/(1-pass$qpq)*1000.0
    
    # Get Cslow from long nutrient cycling solution
    omegas <- aequiln$af*pass$omegafs + aequiln$ar*pass$omegars
    CslowLong_no_priming <- omegas*VLongN_on$equilNPP/pass$decomp_s/(1-pass$qsq)*1000.0
    
    # Calculate nutrient release from wody pool
    NrelwoodVLong <- aequiln$aw*aequiln$nw*VLongN_on$equilNPP*1000.0
    
    # Calculate long term nutrient constraint
    NCHUGH_on <- NConsLong(nfseq_on, a_vec, CpassVLong,
                                            NinL = Nin)
    
    # calculate N and C gaps for priming to occur
    c_into_active <- VLongN_on$equilNPP * aequiln$ar * aequiln$ariz * rhizo_cue * 1000.0
    n_into_active <- c_into_active * aequiln$nr
    n_active_gap <- c_into_active * nca - n_into_active
    
    # adjust decomposition of slow pool to close the N gap
    new_kdec <- pass$decomp_s * (1 + km) * pmax(c_into_active/(c_into_active + km), 0.3)
    
    # Calculate C slow based on exudation and new decomposition values
    CslowLong <- omegas*VLongN_on$equilNPP/new_kdec/(1-pass$qsq)*1000.0
    
    # Calculate medium term nutrient constraint
    NCMEDIUM_on <- NConsMedium_priming(df=nfseq_on, 
                                       a=a_vec, 
                                       Cpass=CpassVLong, 
                                       Cslow=CslowLong, 
                                       NinL = Nin+NrelwoodVLong)
    
    # Solve longterm equilibrium
    equil_long_350_on <- solveLong_full_cn_medium(CO2=CO2_1, Cpass=CpassVLong, NinL = Nin)#+NrelwoodVLong)
    equil_long_700_on <- solveLong_full_cn_medium(CO2=CO2_2, Cpass=CpassVLong, NinL = Nin)#+NrelwoodVLong)
    
    # Solve medium equilibrium
    equil_medium_350_on <- solveMedium_priming(CO2=CO2_1, Cpass=CpassVLong, Cslow=CslowLong, NinL = Nin+NrelwoodVLong)
    equil_medium_700_on <- solveMedium_priming(CO2=CO2_2, Cpass=CpassVLong, Cslow=CslowLong, NinL = Nin+NrelwoodVLong)
    
    # get the point instantaneous NPP response to doubling of CO2
    df700 <- as.data.frame(cbind(round(nfseq_on,3), PC700_on))
    inst700_on <- inst_NPP(VLongN_on$equilnf, df700)
    
    ## locate the intersect between VL nutrient constraint and CO2 = 700
    VLong700_on <- solveVLong_exudation_medium(CO2_2)
    
    out.tab[8,"NPP_350"] <- VLongN_on$equilNPP
    out.tab[8,"NPP_700"] <- VLong700_on$equilNPP
    out.tab[8,"I"] <- (inst700_on$equilNPP - VLongN_on$equilNPP) / VLongN_on$equilNPP
    out.tab[8,"M"] <- (equil_medium_700_on$equilNPP - VLongN_on$equilNPP) / VLongN_on$equilNPP
    out.tab[8,"L"] <- (equil_long_700_on$equilNPP - VLongN_on$equilNPP) / VLongN_on$equilNPP
    out.tab[8,"VL"] <- ( VLong700_on$equilNPP - VLongN_on$equilNPP) / VLongN_on$equilNPP
    
    
    
    
    
    ######## Save the output table
    write.csv(out.tab, paste0(destDir, "/Table1.csv"), row.names=F)
}

######################## Script ###################################
compute_Table_1(destDir="Tables")

