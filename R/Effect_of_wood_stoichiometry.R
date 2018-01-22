

wood_stoichiometry_effect <- function() {
    #### Perform CNP only analysis of VL pools
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
    
    # Calculate pf based on nf of long-term nutrient exchange
    pfseqL <- inferpfL(nfseq, a_nf, PinL = Pin,#+PrelwoodVLong,
                       NinL = Nin,#+NrelwoodVLong,
                       Cpass=CpassVLong)
    
    # Calculate long term nutrieng constraint
    NCLONG_vary <- Long_constraint_N(nfseq, a_nf, CpassVLong,
                                     NinL = Nin)#+NrelwoodVLong)
    
    # Find long term equilibrium point
    Long_equil_vary_350 <- solveLong_full_cnp(CO2=CO2_1, Cpass=CpassVLong, NinL = Nin,#+NrelwoodVLong, 
                                          PinL=Pin)#+PrelwoodVLong)
    Long_equil_vary_700 <- solveLong_full_cnp(CO2=CO2_2, Cpass=CpassVLong, NinL = Nin,#+NrelwoodVLong, 
                                          PinL=Pin)#+PrelwoodVLong)
    
    omegas <- aequiln$af*pass$omegafs + aequiln$ar*pass$omegars
    CslowLong <- omegas*Long_equil_vary_350$equilNPP/pass$decomp_s/(1-pass$qsq)*1000.0
    
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
    
    Medium_equil_vary_350 <- solveMedium_full_cnp(CO2=CO2_1, Cpass=CpassVLong, Cslow=CslowLong, NinL = Nin+NrelwoodVLong,
                                                  PinL=Pin+PrelwoodVLong)
    Medium_equil_vary_700 <- solveMedium_full_cnp(CO2=CO2_2, Cpass=CpassVLong, Cslow=CslowLong, NinL = Nin+NrelwoodVLong,
                                                  PinL=Pin+PrelwoodVLong)
    
    #### Perform CNP only analysis of VL pools
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
    
    # Calculate pf based on nf of long-term nutrient exchange
    pfseqL <- inferpfL(nfseq, a_nf, PinL = Pin,#+PrelwoodVLong,
                       NinL = Nin,#+NrelwoodVLong,
                       Cpass=CpassVLong)
    
    # Calculate long term nutrieng constraint
    NCLONG_fix <- Long_constraint_N(nfseq, a_nf, CpassVLong,
                                     NinL = Nin)#+NrelwoodVLong)
    
    # Find long term equilibrium point
    Long_equil_fix_350 <- solveLong_full_cnp(CO2=CO2_1, Cpass=CpassVLong, NinL = Nin,#+NrelwoodVLong, 
                                          PinL=Pin)#+PrelwoodVLong)
    Long_equil_fix_700 <- solveLong_full_cnp(CO2=CO2_2, Cpass=CpassVLong, NinL = Nin,#+NrelwoodVLong, 
                                         PinL=Pin)#+PrelwoodVLong)
    
    omegas <- aequiln$af*pass$omegafs + aequiln$ar*pass$omegars
    CslowLong <- omegas*Long_equil_fix_350$equilNPP/pass$decomp_s/(1-pass$qsq)*1000.0
    
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
    
    Medium_equil_fix_350 <- solveMedium_full_cnp(CO2=CO2_1, Cpass=CpassVLong, Cslow=CslowLong, NinL = Nin+NrelwoodVLong,
                                                  PinL=Pin+PrelwoodVLong)
    Medium_equil_fix_700 <- solveMedium_full_cnp(CO2=CO2_2, Cpass=CpassVLong, Cslow=CslowLong, NinL = Nin+NrelwoodVLong,
                                                  PinL=Pin+PrelwoodVLong)
    
    
    # CO2 effect
    co2_effect_vary <- (VLong_equil_vary_700$equilNPP - VLong_equil_vary_350$equilNPP) / VLong_equil_vary_350$equilNPP * 100
    co2_effect_fix <- (VLong_equil_fix_700$equilNPP - VLong_equil_fix_350$equilNPP) / VLong_equil_fix_350$equilNPP * 100

    #### Plotting
    tiff("Plots/Effect_of_wood_stoichiometry.tiff",
         width = 10, height = 5, units = "in", res = 300)
    par(mfrow=c(1,2), mar=c(5.1,6.1,2.1,2.1))
    
    # shoot nc vs. NPP at variable wood
    plot(nfseq, photo_350_vary, xlim=c(0.005, 0.03),
         ylim=c(0.5, 3), 
         type = "l", xlab = "Shoot N:C ratio", 
         ylab = expression(paste("Production [kg C ", m^-2, " ", yr^-1, "]")),
         col="cyan", lwd = 3, cex.lab=1.5)
    points(nfseq, vlong_vary$NPP_N, type="l", col="tomato", lwd = 3)
    points(VLong_equil_vary_350$equilnf, VLong_equil_vary_350$equilNPP, type="p", pch = 19, col = "blue", cex = 2)
    points(nfseq, NCLONG_vary$NPP, type='l',col="violet", lwd = 3)
    points(nfseq, photo_700_vary, col="green", type="l", lwd = 3)
    points(Long_equil_vary_700$equilnf, Long_equil_vary_700$equilNPP, type="p", col="orange", pch = 19, cex = 2)
    points(VLong_equil_vary_700$equilnf, VLong_equil_vary_700$equilNPP,type="p", col="red", pch = 19, cex = 2)
    points(Medium_equil_vary_700$equilnf, Medium_equil_vary_700$equilNPP, type="p", col="purple", pch = 19, cex = 2)
    
    text(x=0.006, y=2.9, "(a)", cex = 2)
    
    # shoot nc vs. NPP at fixed wood
    plot(nfseq, photo_350_fix, xlim=c(0.005, 0.03),
         ylim=c(0.5, 3), 
         type = "l", xlab = "Shoot N:C ratio", 
         ylab = expression(paste("Production [kg C ", m^-2, " ", yr^-1, "]")),
         col="cyan", lwd = 3, cex.lab=1.5)
    points(nfseq, vlong_fix$NPP_N, type="l", col="tomato", lwd = 3)
    points(VLong_equil_fix_350$equilnf, VLong_equil_fix_350$equilNPP, type="p", pch = 19, col = "blue", cex = 2)
    points(nfseq, NCLONG_fix$NPP, type='l',col="violet", lwd = 3)
    points(nfseq, photo_700_fix, col="green", type="l", lwd = 3)
    points(Long_equil_fix_700$equilnf, Long_equil_fix_700$equilNPP, type="p", col="orange", pch = 19, cex = 2)
    points(VLong_equil_fix_700$equilnf, VLong_equil_fix_700$equilNPP,type="p", col="red", pch = 19, cex = 2)
    points(Medium_equil_fix_700$equilnf, Medium_equil_fix_700$equilNPP, type="p", col="purple", pch = 19, cex = 2)
    
    text(x=0.006, y=2.9, "(b)", cex = 2)
    

    legend("topright", c("P350", "P700", "VL", "L",
                         "A", "B", "C", "D"),
           col=c("cyan","green", "tomato", "violet","blue", "purple", "red", "orange"), 
           lwd=c(2,2,2,2,NA,NA,NA,NA), pch=c(NA,NA,NA, NA,19,19,19,19), cex = 1.2, 
           bg = adjustcolor("grey", 0.8), ncol=2)
    
    dev.off()
    
}


#### Program
wood_stoichiometry_effect()