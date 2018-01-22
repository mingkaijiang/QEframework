
#### Functions to generate Table on wood stoichiometric flexibility
#### Purpose:
#### To compute pool sizes of slow SOM, passive SOM and wood biomass
#### comparing fixed vs. variable wood stoichiometry
#### to demonstrate the effect of slow SOM on CO2 fertilization effect
####
################################################################################

######################## Functions ###################################
compute_stoichiometric_flexibility_table <- function(destDir) {

    
    ######## Make the output table
    out.tab <- matrix(ncol=10, nrow = 2)
    colnames(out.tab) <- c("model", "Cpass", "Npass", "Ppass",
                           "Cslow", "Nslow", "Pslow",
                           "Cwood", "Nwood", "Pwood")
    out.tab <- as.data.frame(out.tab)
    out.tab[1,"model"] <- "Variable"
    out.tab[2,"model"] <- "Fixed"
    
    
    #### Analytical run 1
    source("Parameters/Analytical_Run1_Parameters.R")
    
    # create a range of nc for shoot to initiate
    nfseq <- round(seq(0.001, 0.1, by = 0.001),5)
    a_nf <- as.data.frame(allocn(nfseq))
    
    # using very long term relationship to calculate pf from nf
    pfseq <- inferpfVL(nfseq, a_nf)
    a_pf <- as.data.frame(allocp(pfseq))
    
    # calculate photosynthetic constraint at CO2 = 350
    photo_350_vary <- photo_constraint_full_cnp(nfseq, pfseq, a_nf, a_pf, CO2_1)
    photo_700_vary <- photo_constraint_full_cnp(nfseq, pfseq, a_nf, a_pf, CO2_2)
    
    ### calculate very long term NC and PC constraint on NPP, respectively
    vlong_vary <- VLong_constraint_N(nf=nfseq, nfdf=a_nf)
    
    ### finding the equilibrium point between photosynthesis and very long term nutrient constraints
    VLong_equil_vary_350 <- solveVLong_full_cnp(CO2=CO2_1)
    VLong_equil_vary_700 <- solveVLong_full_cnp(CO2=CO2_2)
    
    ### Get Cpassive from very-long nutrient cycling solution
    aequiln <- allocn(VLong_equil_vary_350$equilnf)
    aequilp <- allocp(VLong_equil_vary_350$equilpf)
    
    pass <- slow_pool(df=VLong_equil_vary_350$equilnf, a=aequiln)
    omegap <- aequiln$af*pass$omegafp + aequiln$ar*pass$omegarp
    CpassVLong <- omegap*VLong_equil_vary_350$equilNPP/pass$decomp_p/(1-pass$qpq)*1000.0
    
    ### Calculate nutrient release from recalcitrant pools
    PrelwoodVLong <- aequilp$aw*aequilp$pw*VLong_equil_vary_350$equilNPP*1000.0
    NrelwoodVLong <- aequiln$aw*aequiln$nw*VLong_equil_vary_350$equilNPP*1000.0
    
    # Find long term equilibrium point
    Long_equil_vary_350 <- solveLong_full_cnp(CO2=CO2_1, Cpass=CpassVLong, NinL = Nin,#+NrelwoodVLong, 
                                              PinL=Pin)#+PrelwoodVLong)
    
    omegas <- aequiln$af*pass$omegafs + aequiln$ar*pass$omegars
    CslowLong <- omegas*Long_equil_vary_350$equilNPP/pass$decomp_s/(1-pass$qsq)*1000.0
    
    out.tab[1,"Cpass"] <- CpassVLong
    out.tab[1,"Cslow"] <- CslowLong
    out.tab[1,"Cwood"] <- aequilp$aw*VLong_equil_vary_350$equilNPP*1000.0
    
    out.tab[1,"Npass"] <- CpassVLong * ncp
    out.tab[1,"Nslow"] <- CslowLong * ncs
    out.tab[1,"Nwood"] <- NrelwoodVLong
    
    out.tab[1,"Ppass"] <- CpassVLong * pcp
    out.tab[1,"Pslow"] <- CslowLong * pcs
    out.tab[1,"Pwood"] <- PrelwoodVLong
    
    
    #### analytical run 3
    source("Parameters/Analytical_Run3_Parameters.R")
    
    # create a range of nc for shoot to initiate
    nfseq <- round(seq(0.001, 0.1, by = 0.001),5)
    a_nf <- as.data.frame(allocn(nfseq))
    
    # using very long term relationship to calculate pf from nf
    pfseq <- inferpfVL(nfseq, a_nf)
    a_pf <- as.data.frame(allocp(pfseq))
    
    # calculate photosynthetic constraint at CO2 = 350
    photo_350_fix <- photo_constraint_full_cnp(nfseq, pfseq, a_nf, a_pf, CO2_1)
    photo_700_fix <- photo_constraint_full_cnp(nfseq, pfseq, a_nf, a_pf, CO2_2)
    
    ### calculate very long term NC and PC constraint on NPP, respectively
    vlong_fix <- VLong_constraint_N(nf=nfseq, nfdf=a_nf)
    
    ### finding the equilibrium point between photosynthesis and very long term nutrient constraints
    VLong_equil_fix_350 <- solveVLong_full_cnp(CO2=CO2_1)
    VLong_equil_fix_700 <- solveVLong_full_cnp(CO2=CO2_2)
    
    ### Get Cpassive from very-long nutrient cycling solution
    aequiln <- allocn(VLong_equil_fix_350$equilnf)
    aequilp <- allocp(VLong_equil_fix_350$equilpf)
    pass <- slow_pool(df=VLong_equil_fix_350$equilnf, a=aequiln)
    omegap <- aequiln$af*pass$omegafp + aequiln$ar*pass$omegarp
    CpassVLong <- omegap*VLong_equil_fix_350$equilNPP/pass$decomp_p/(1-pass$qpq)*1000.0
    
    ### Calculate nutrient release from recalcitrant pools
    PrelwoodVLong <- aequilp$aw*aequilp$pw*VLong_equil_fix_350$equilNPP*1000.0
    NrelwoodVLong <- aequiln$aw*aequiln$nw*VLong_equil_fix_350$equilNPP*1000.0
    
    # Find long term equilibrium point
    Long_equil_fix_350 <- solveLong_full_cnp(CO2=CO2_1, Cpass=CpassVLong, NinL = Nin,#+NrelwoodVLong, 
                                             PinL=Pin)#+PrelwoodVLong)
    
    omegas <- aequiln$af*pass$omegafs + aequiln$ar*pass$omegars
    CslowLong <- omegas*Long_equil_fix_350$equilNPP/pass$decomp_s/(1-pass$qsq)*1000.0
    

    out.tab[2,"Cpass"] <- CpassVLong
    out.tab[2,"Cslow"] <- CslowLong
    out.tab[2,"Cwood"] <- aequilp$aw*VLong_equil_fix_350$equilNPP*1000.0
    
    out.tab[2,"Npass"] <- CpassVLong * ncp
    out.tab[2,"Nslow"] <- CslowLong * ncs
    out.tab[2,"Nwood"] <- NrelwoodVLong
    
    out.tab[2,"Ppass"] <- CpassVLong * pcp
    out.tab[2,"Pslow"] <- CslowLong * pcs
    out.tab[2,"Pwood"] <- PrelwoodVLong

    
    ######## Save the output table
    write.csv(out.tab, paste0(destDir, "/Stoichiometric_flexibility_table.csv"), row.names=F)
    
}


######################## Script ###################################
compute_stoichiometric_flexibility_table(destDir="Tables")