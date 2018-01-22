
#### Analytical script Run 1
####
#### Assumptions:
#### 1. baseline model
#### 2. Variable wood NC
#### 3. baseline N cycle
####
################################################################################
#### Functions
Perform_Analytical_Run1 <- function(f.flag = 1, cDF, eDF) {
    #### Function to perform analytical run 1 simulations
    #### eDF: stores equilibrium points
    #### cDF: stores constraint points (curves)
    #### f.flag: = 1 simply plot analytical solution file
    #### f.flag: = 2 return cDF
    #### f.flag: = 3 return eDF

    ######### Main program
    source("Parameters/Analytical_Run1_Parameters.R")
    
    ### create a range of nc for shoot to initiate
    nfseq <- round(seq(0.001, 0.1, by = 0.001),5)
    a_nf <- as.data.frame(allocn(nfseq))
    
    ### calculate photosynthetic constraint at CO2 = 350
    Photo350 <- photo_constraint_full_cn(nfseq, a_nf, CO2_1)

    ### calculate very long term NC and PC constraint on NPP, respectively
    NCVLONG <- VLong_constraint_N(nf=nfseq, nfdf=a_nf)
    
    ### finding the equilibrium point between photosynthesis and very long term nutrient constraints
    VLong_equil <- solveVLong_full_cn(CO2=CO2_1)
    
    ### Get Cpassive from very-long nutrient cycling solution
    aequiln <- allocn(VLong_equil$equilnf)
    pass <- slow_pool(df=VLong_equil$equilnf, a=aequiln)
    omegap <- aequiln$af*pass$omegafp + aequiln$ar*pass$omegarp
    CpassVLong <- omegap*VLong_equil$equilNPP/pass$decomp_p/(1-pass$qpq)*1000.0
    
    # Calculate long term nutrient constraint
    NCLONG <- Long_constraint_N(nfseq, a_nf, CpassVLong,
                                NinL = Nin)#+NrelwoodVLong)
    
    # Find long term equilibrium point
    Long_equil <- solveLong_full_cn(CO2=CO2_1, Cpass=CpassVLong, NinL = Nin)#+NrelwoodVLong)
    
    # Get Cslow from long nutrient cycling solution
    omegas <- aequiln$af*pass$omegafs + aequiln$ar*pass$omegars
    CslowLong <- omegas*Long_equil$equilNPP/pass$decomp_s/(1-pass$qsq)*1000.0
    
    ### Calculate nutrient release from slow woody pool
    NrelwoodVLong <- aequiln$aw*aequiln$nw*VLong_equil$equilNPP*1000.0
    
    # Calculate medium term nutrient constraint
    NCMEDIUM <- NConsMedium(df=nfseq, 
                            a=a_nf, 
                            Cpass=CpassVLong, 
                            Cslow=CslowLong, 
                            NinL = Nin+NrelwoodVLong)
    
    Medium_equil_350 <- solveMedium(CO2=CO2_1, 
                                    Cpass=CpassVLong, 
                                    Cslow=CslowLong, 
                                    NinL = Nin+NrelwoodVLong)
    

    out350DF <- data.frame(nfseq, Photo350, NCVLONG, NCLONG)
    colnames(out350DF) <- c("nc", "NPP_350", "NPP_VL",
                            "nleach_VL", "NPP_350_L", "nwood_L", "nburial_L",
                            "nleach_L", "aw")
    equil350DF <- data.frame(VLong_equil, Long_equil)
    colnames(equil350DF) <- c("nc_VL", "NPP_VL", 
                              "nc_L", "NPP_L")
    
    # store constraint and equil DF onto their respective output df
    # cDF[cDF$Run == 1 & cDF$CO2 == 350, 3:13] <- out350DF
    # eDF[eDF$Run == 1 & eDF$CO2 == 350, 3:8] <- equil350DF
    
    ##### CO2 = 700
    # N:C ratio
    nfseq <- round(seq(0.001, 0.1, by = 0.001),5)
    a_nf <- as.data.frame(allocn(nfseq))
    
    # calculate NC vs. NPP at CO2 = 350 respectively
    Photo700 <- photo_constraint_full_cn(nfseq, a_nf, CO2_2)
    
    ### calculate very long term NC and PC constraint on NPP, respectively
    NCVLONG <- VLong_constraint_N(nf=nfseq, nfdf=a_nf)
    
    ### finding the equilibrium point between photosynthesis and very long term nutrient constraints
    VLong_equil <- solveVLong_full_cn(CO2=CO2_2)
    
    # Find long term equilibrium point
    Long_equil <- solveLong_full_cn(CO2=CO2_2, Cpass=CpassVLong, NinL = Nin)
    
    # Find medium term equilibrium point
    Medium_equil_350 <- solveMedium(CO2_1, Cpass = CpassVLong, Cslow = CslowLong, 
                                             NinL=Nin+NrelwoodVLong)
    Medium_equil_700 <- solveMedium(CO2_2, Cpass = CpassVLong, Cslow = CslowLong, 
                                             NinL=Nin+NrelwoodVLong)
    
    out700DF <- data.frame(nfseq, Photo700, NCVLONG, NCLONG)
    colnames(out700DF) <- c("nc", "NPP_700", "NPP_VL",
                            "nleach_VL", "NPP_700_L", "nwood_L", "nburial_L",
                            "nleach_L", "aw")
    
    equil700DF <- data.frame(VLong_equil, Long_equil)
    colnames(equil700DF) <- c("nc_VL", "NPP_VL", 
                              "nc_L", "NPP_L")
    
    
    # store constraint and equil DF onto their respective output df
    # cDF[cDF$Run == 1 & cDF$CO2 == 700, 3:13] <- out700DF
    # eDF[eDF$Run == 1 & eDF$CO2 == 700, 3:8] <- equil700DF
    
    # get the point instantaneous NPP response to doubling of CO2
    df700 <- as.data.frame(cbind(round(nfseq,3), Photo700))
    inst700 <- inst_NPP(equil350DF$nc_VL, df700)
    
    # eDF[eDF$Run == 1 & eDF$CO2 == 350, 9] <- inst700$equilNPP
    # eDF[eDF$Run == 1 & eDF$CO2 == 700, 9] <- inst700$equilNPP
    
    if (f.flag == 1) {
        
        ### plot 2-d plots of nf vs. npp and nf vs. pf
        tiff("Plots/Analytical_Run1_2d.tiff",
             width = 10, height = 5, units = "in", res = 300)
        par(mfrow=c(1,2), mar=c(5.1,6.1,2.1,2.1))
        
        # shoot nc vs. NPP
        plot(out350DF$nc, out350DF$NPP_350, xlim=c(0.01, 0.05),
              ylim=c(0.5, 3), 
             type = "l", xlab = "Shoot N:C ratio", 
             ylab = expression(paste("Production [kg C ", m^-2, " ", yr^-1, "]")),
             col="cyan", lwd = 3, cex.lab=1.5)
        points(out350DF$nc, out350DF$NPP_VL, type="l", col="tomato", lwd = 3)
        points(equil350DF$nc_VL, equil350DF$NPP_VL, type="p", pch = 19, col = "blue", cex = 2)
        points(out350DF$nc, out350DF$NPP_350_L, type='l',col="violet", lwd = 3)
        
        points(nfseq, NCMEDIUM$NPP, type="l", col="darkred", lwd = 3)
        
        points(out700DF$nc, out700DF$NPP_700, col="green", type="l", lwd = 3)
        points(equil350DF$nc_VL, inst700$equilNPP, type="p", col = "darkgreen", pch=19, cex = 2)
        points(equil700DF$nc_VL, equil700DF$NPP_VL, type="p", col="orange", pch = 19, cex = 2)
        points(equil700DF$nc_L, equil700DF$NPP_L,type="p", col="red", pch = 19, cex = 2)
        points(Medium_equil_700$equilnf, Medium_equil_700$equilNPP, type="p", col="purple", pch = 19, cex = 2)
        text(x=0.045, y=2.9, "(a)", cex = 2)
        
        plot(0,0)
        legend("bottomright", c("P350", "P700", "VL", "L", "M",
                            "A", "B", "C", "D", "E"),
               col=c("cyan","green", "tomato", "violet","darkred","blue", "darkgreen","purple","red", "orange"), 
               lwd=c(2,2,2,2,2,NA,NA,NA,NA,NA), pch=c(NA,NA,NA,NA,NA,19,19,19,19,19), cex = 1.0, 
               bg = adjustcolor("grey", 0.8), ncol=2)
        
        dev.off()
        
    } else if (f.flag == 2) {
        return(cDF)
    } else if (f.flag == 3) {
        return(eDF)
    }
    
}
