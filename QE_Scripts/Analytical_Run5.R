
#### Analytical script to match GDAY Run 5 settings
####
#### Assumptions:
#### Same as Run 1, but with
#### 1. simplified LUE model depends on CNP

################################################################################

#### Functions
Perform_Analytical_Run5 <- function(f.flag = 1, cDF, eDF) {
    #### Function to perform analytical run 5 simulations
    #### eDF: stores equilibrium points
    #### cDF: stores constraint points (curves)
    #### f.flag: = 1 simply plot analytical solution file
    #### f.flag: = 2 return cDF
    #### f.flag: = 3 return eDF

    ######### Main program
    source("Parameters/Analytical_Run5_Parameters.R")
    
    # create a range of nc for shoot to initiate
    nfseq <- round(seq(0.01, 0.1, by = 0.001),5)
    a_nf <- as.data.frame(allocn(nfseq))
    
    # using very long term relationship to calculate pf from nf
    pfseq <- inferpfVL(nfseq, a_nf)
    a_pf <- as.data.frame(allocp(pfseq))
    
    # calculate photosynthetic constraint at CO2 = 350
    Photo350 <- photo_constraint_simple_cnp(nfseq, pfseq, a_nf, a_pf, CO2_1)
    
    ### calculate very long term NC and PC constraint on NPP, respectively
    NCVLONG <- VLong_constraint_N(nf=nfseq, nfdf=a_nf)
    
    ### NPP derived from PCVLONG should match NPP from NCVLONG
    PCVLONG <- VLong_constraint_P(pf=pfseq, pfdf=a_pf)
    
    ### finding the equilibrium point between photosynthesis and very long term nutrient constraints
    VLong_equil <- solveVLong_simple_cnp(CO2=CO2_1)
    
    ### Get Cpassive from very-long nutrient cycling solution
    aequiln <- allocn(VLong_equil$equilnf)
    aequilp <- allocp(VLong_equil$equilpf)
    pass <- passive(df=VLong_equil$equilnf, a=aequiln)
    omega <- aequiln$af*pass$omegaf + aequiln$ar*pass$omegar
    CpassVLong <- omega*VLong_equil$equilNPP/pass$decomp/(1-pass$qq)*1000.0
    
    ### Calculate nutrient release from recalcitrant pools
    PrelwoodVLong <- aequilp$aw*aequilp$pw*VLong_equil$equilNPP*1000.0
    NrelwoodVLong <- aequiln$aw*aequiln$nw*VLong_equil$equilNPP*1000.0
    
    # Calculate pf based on nf of long-term nutrient exchange
    pfseqL <- inferpfL(nfseq, a_nf, PinL = Pin,#+PrelwoodVLong,
                       NinL = Nin,#+NrelwoodVLong,
                       Cpass=CpassVLong)
    
    # Calculate long term nutrieng constraint
    NCLONG <- Long_constraint_N(nfseq, a_nf, CpassVLong,
                                NinL = Nin)#+NrelwoodVLong)
    
    PCLONG <- Long_constraint_P(nfseq, pfseqL, allocp(pfseqL),
                                CpassVLong, PinL=Pin)#+PrelwoodVLong)
    
    # Find long term equilibrium point
    Long_equil <- solveLong_simple_cnp(CO2=CO2_1, Cpass=CpassVLong, NinL = Nin,#+NrelwoodVLong, 
                                     PinL=Pin)#+PrelwoodVLong)
    
    
    out350DF <- data.frame(nfseq, pfseq, pfseqL, Photo350, NCVLONG, NCLONG)
    colnames(out350DF) <- c("nc", "pc_VL", "pc_350_L", "NPP_350", "NPP_VL",
                            "nleach_VL", "NPP_350_L", "nwood_L", "nburial_L",
                            "nleach_L", "aw")
    equil350DF <- data.frame(VLong_equil, Long_equil)
    colnames(equil350DF) <- c("nc_VL", "pc_VL", "NPP_VL", 
                              "nc_L","pc_L", "NPP_L")
    
    # store constraint and equil DF onto their respective output df
    #cDF[cDF$Run == 5 & cDF$CO2 == 350, 3:13] <- out350DF
    #eDF[eDF$Run == 5 & eDF$CO2 == 350, 3:8] <- equil350DF
    
    ##### CO2 = 700
    
    # N:C and P:C ratio
    nfseq <- round(seq(0.01, 0.1, by = 0.001),5)
    a_nf <- as.data.frame(allocn(nfseq))
    
    # using very long term relationship to calculate pf from nf
    pfseq <- inferpfVL(nfseq, a_nf)
    a_pf <- as.data.frame(allocp(pfseq))
    
    # calculate NC vs. NPP at CO2 = 700 respectively
    Photo700 <- photo_constraint_simple_cnp(nfseq, pfseq, a_nf, a_pf, CO2_2)
    
    ### calculate very long term NC and PC constraint on NPP, respectively
    NCVLONG <- VLong_constraint_N(nf=nfseq, nfdf=a_nf)
    
    ### NPP derived from PCVLONG should match NPP from NCVLONG
    PCVLONG <- VLong_constraint_P(pf=pfseq, pfdf=a_pf)
    
    ### finding the equilibrium point between photosynthesis and very long term nutrient constraints
    VLong_equil <- solveVLong_simple_cnp(CO2=CO2_2)
    
    # Find long term equilibrium point
    Long_equil <- solveLong_simple_cnp(CO2=CO2_2, Cpass=CpassVLong, NinL = Nin,#+NrelwoodVLong, 
                                     PinL=Pin)#+PrelwoodVLong)
    
    out700DF <- data.frame(nfseq, pfseq, pfseqL, Photo700, NCVLONG, NCLONG)
    colnames(out700DF) <- c("nc", "pc_VL", "pc_700_L", "NPP_700", "NPP_VL",
                            "nleach_VL", "NPP_700_L", "nwood_L", "nburial_L",
                            "nleach_L", "aw")
    
    equil700DF <- data.frame(VLong_equil, Long_equil)
    colnames(equil700DF) <- c("nc_VL", "pc_VL", "NPP_VL", 
                              "nc_L","pc_L", "NPP_L")
    
    
    # store constraint and equil DF onto their respective output df
    #cDF[cDF$Run == 5 & cDF$CO2 == 700, 3:13] <- out700DF
    #eDF[eDF$Run == 5 & eDF$CO2 == 700, 3:8] <- equil700DF
    
    # get the point instantaneous NPP response to doubling of CO2
    df700 <- as.data.frame(cbind(round(nfseq,3), Photo700))
    inst700 <- inst_NPP(equil350DF$nc_VL, df700)
    
    #eDF[eDF$Run == 5 & eDF$CO2 == 350, 9] <- inst700$equilNPP
    #eDF[eDF$Run == 5 & eDF$CO2 == 700, 9] <- inst700$equilNPP
    
    if (f.flag == 1) {
        
        
        #### Library
        require(scatterplot3d)
        
        ######### Plotting
        
        tiff("Plots/Analytical_Run5.tiff",
             width = 8, height = 7, units = "in", res = 300)
        par(mar=c(5.1,5.1,2.1,2.1))
        
        
        # NPP constraint by CO2 = 350
        s3d <- scatterplot3d(out350DF$nc, out350DF$pc_VL, out350DF$NPP_350, xlim=c(0.0, 0.1),
                             ylim = c(0.0, 0.02), zlim=c(0, 5), 
                             type = "l", xlab = "Shoot N:C ratio", ylab = "Shoot P:C ratio", 
                             zlab = expression(paste("Production [kg C ", m^-2, " ", yr^-1, "]")),
                             color="cyan", lwd = 3, angle=24)
        
        # NPP constraint by very long term nutrient availability
        s3d$points3d(out350DF$nc, out350DF$pc_VL, out350DF$NPP_VL, type="l", col="tomato", lwd = 3)
        
        # equilibrated NPP for very long term nutrient and CO2 = 350
        s3d$points3d(equil350DF$nc_VL, equil350DF$pc_VL, equil350DF$NPP_VL,
                     type="h", pch = 19, col = "blue")
        
        # NPP constraint by long term nutrient availability
        s3d$points3d(out350DF$nc, out350DF$pc_VL, out350DF$NPP_350_L, type='l',col="violet", lwd = 3)
        #s3d$points3d(out700DF$nc, out700DF$pc_700_L, out700DF$NPP_700_L, type='l',col="grey", lwd = 3)
        
        
        # equilibrated NPP for long term nutrient and CO2 = 350
        #s3d$points3d(equil350DF$nc_L, equil350DF$pc_L, equil350DF$NPP_L,
        #             type="h", col="lightblue", pch = 19)
        
        # NPP constraint by CO2 = 700
        s3d$points3d(out700DF$nc, out700DF$pc_VL, out700DF$NPP_700, col="green", type="l", lwd = 3)
        
        s3d$points3d(equil350DF$nc_VL, equil350DF$pc_VL, 
                     inst700$equilNPP, type="h", col = "darkgreen", pch=19)
        
        # equilibrated NPP for very long term nutrient and CO2 = 700
        s3d$points3d(equil700DF$nc_VL, equil700DF$pc_VL, equil700DF$NPP_VL, 
                     type="h", col="orange", pch = 19)
        
        # equilibrated NPP for long term nutrient and CO2 = 700
        s3d$points3d(equil700DF$nc_L, equil700DF$pc_VL, equil700DF$NPP_L,
                     type="h", col="red", pch = 19)
        
        
        legend("topleft", c(expression(paste("Photo constraint at ", CO[2]," = 350 ppm")), 
                            expression(paste("Photo constraint at ", CO[2]," = 700 ppm")), 
                            "VL nutrient constraint", "L nutrient constraint",
                            "A", "B", "C", "D"),
               col=c("cyan","green", "tomato", "violet","blue", "darkgreen","red", "orange"), 
               lwd=c(2,2,2,2,NA,NA,NA,NA), pch=c(NA,NA,NA,NA,19,19,19,19), cex = 1.0, 
               bg = adjustcolor("grey", 0.8))
        
        dev.off()
        
        ### plot 2-d plots of nf vs. npp and nf vs. pf
        tiff("Plots/Analytical_Run5_2d.tiff",
             width = 10, height = 5, units = "in", res = 300)
        par(mfrow=c(1,2), mar=c(5.1,6.1,2.1,2.1))
        
        # shoot nc vs. NPP
        plot(out350DF$nc, out350DF$NPP_350, xlim=c(0.0, 0.1),
             ylim=c(0, 3), 
             type = "l", xlab = "Shoot N:C ratio", 
             ylab = expression(paste("Production [kg C ", m^-2, " ", yr^-1, "]")),
             col="cyan", lwd = 3,cex.lab=1.5)
        points(out350DF$nc, out350DF$NPP_VL, type="l", col="tomato", lwd = 3)
        points(equil350DF$nc_VL, equil350DF$NPP_VL, type="p", pch = 19, col = "blue",cex=2)
        points(out350DF$nc, out350DF$NPP_350_L, type='l',col="violet", lwd = 3)
        points(out700DF$nc, out700DF$NPP_700, col="green", type="l", lwd = 3)
        points(equil350DF$nc_VL, inst700$equilNPP, type="p", col = "darkgreen", pch=19,cex=2)
        points(equil700DF$nc_VL, equil700DF$NPP_VL, type="p", col="orange", pch = 19,cex=2)
        points(equil700DF$nc_L, equil700DF$NPP_L,type="p", col="red", pch = 19,cex=2)

        
        # shoot nc vs. shoot pc
        plot(out350DF$nc, out350DF$pc_VL, xlim=c(0.0, 0.1),
             ylim=c(0, 0.02), 
             type = "l", xlab = "Shoot N:C ratio", 
             ylab = "Shoot P:C ratio",
             col="cyan", lwd = 3,cex.lab=1.5)
        points(out350DF$nc, out350DF$pc_VL, type="l", col="tomato", lwd = 3)
        
        points(equil350DF$nc_VL, equil350DF$pc_VL, type="p", pch = 19, col = "blue",cex=2)
        
        points(out350DF$nc, out350DF$pc_VL, type='l',col="violet", lwd = 3)
        
        points(out700DF$nc, out700DF$pc_VL, col="green", type="l", lwd = 3)
        
        points(equil350DF$nc_VL, equil350DF$pc_VL, type="p", col = "darkgreen", pch=19,cex=2)
        
        points(equil700DF$nc_VL, equil700DF$pc_VL, type="p", col="orange", pch = 19,cex=2)
        
        points(equil700DF$nc_L, equil700DF$pc_VL, type="p", col="red", pch = 19,cex=2)
        
        legend("topright", c(expression(paste("Photo constraint at ", CO[2]," = 350 ppm")), 
                            expression(paste("Photo constraint at ", CO[2]," = 700 ppm")), 
                            "VL nutrient constraint", "L nutrient constraint",
                            "A", "B", "C", "D"),
               col=c("cyan","green", "tomato", "violet","blue", "darkgreen","red", "orange"), 
               lwd=c(2,2,2,2,NA,NA,NA,NA), pch=c(NA,NA,NA,NA,19,19,19,19), cex = 0.7, 
               bg = adjustcolor("grey", 0.8))
        
        
        dev.off()
        
    } else if (f.flag == 2) {
        return(cDF)
    } else if (f.flag == 3) {
        return(eDF)
    }
}
