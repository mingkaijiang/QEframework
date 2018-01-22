
#### Analytical script to match GDAY Run 7 settings
####
#### Same as Run 1, except
#### 1. N mineral pool is explicitly simulated
#### 2. Added a N uptake rate coefficient (i.e. adjustable parameter)
#### 
################################################################################


#### Functions
Perform_Analytical_Run7 <- function(f.flag = 1, cDF, eDF) {
    #### Function to perform analytical run 1 simulations
    #### eDF: stores equilibrium points
    #### cDF: stores constraint points (curves)
    #### f.flag: = 1 simply plot analytical solution file
    #### f.flag: = 2 return cDF
    #### f.flag: = 3 return eDF

    ######### Main program
    source("Parameters/Analytical_Run7_Parameters.R")
    
    # create nc and pc for shoot to initiate
    nfseq <- round(seq(0.001, 0.1, by = 0.001),5)
    a_nf <- as.data.frame(allocn(nfseq))
    
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
    NCMEDIUM <- NConsMedium_expl_min(nfseq, a_nf,CpassVLong, CslowLong,
                                     NinL = Nin+NrelwoodVLong)
    
    out350DF <- data.frame(nfseq, NC350, NCVLONG, NCHUGH)
    colnames(out350DF) <- c("nc", "NPP_350", "NPP_VL",
                            "nleach_VL", "NPP_350_L", "nwood_L", "nburial_L",
                            "nleach_L", "aw")
    equil350DF <- data.frame(VLongN, LongN)
    colnames(equil350DF) <- c("nc_VL", "NPP_VL", 
                              "nc_L", "NPP_L")
    
    ##### CO2 = 700
    # N:C and P:C ratio
    nfseq <- round(seq(0.001, 0.1, by = 0.001),5)
    a_nf <- as.data.frame(allocn(nfseq))
    
    # calculate NC vs. NPP at CO2 = 350 respectively
    NC700 <- photo_constraint_full_cn(nfseq, a_nf, CO2_2)
    
    # calculate very long term NC and PC constraint on NPP, respectively
    NCVLONG <- NConsVLong_expl_min(df=nfseq,a=a_nf)
    
    # solve very-long nutrient cycling constraint
    VLongN <- solveVLong_expl_min(CO2_2)
    
    out700DF <- data.frame(nfseq, NC700, NCVLONG, NCHUGH)
    colnames(out700DF) <- c("nc", "NPP_700", "NPP_VL",
                            "nleach_VL", "NPP_700_L", "nwood_L", "nburial_L",
                            "nleach_L", "aw")
    
    # Find equilibrate intersection and plot
    LongN <- solveLong_expl_min(CO2_2, Cpass=CpassVLong, NinL=Nin)#+NrelwoodVLong)
    
    equil700DF <- data.frame(VLongN, LongN)
    colnames(equil700DF) <- c("nc_VL", "NPP_VL", 
                              "nc_L", "NPP_L")
    
    # Find medium term equilibrium point
    Medium_equil_350 <- solveMedium_expl_min(CO2_1, Cpass = CpassVLong, Cslow = CslowLong, 
                                             NinL=Nin+NrelwoodVLong)
    Medium_equil_700 <- solveMedium_expl_min(CO2_2, Cpass = CpassVLong, Cslow = CslowLong, 
                                             NinL=Nin+NrelwoodVLong)
    
    # get the point instantaneous NPP response to doubling of CO2
    df700 <- as.data.frame(cbind(round(nfseq,3), NC700))
    inst700 <- inst_NPP(equil350DF$nc_VL, df700)
    
    if (f.flag == 1) {
        
        ######### Plotting
        
        ### plot 2-d plots of nf vs. npp and nf vs. pf
        tiff("Plots/Analytical_Run7_2d.tiff",
             width = 10, height = 5, units = "in", res = 300)

        # shoot nc vs. NPP
        plot(out350DF$nc, out350DF$NPP_350, xlim=c(0.0, 0.05),
             ylim=c(0, 3), 
             type = "l", xlab = "Shoot N:C ratio", 
             ylab = expression(paste("Production [kg C ", m^-2, " ", yr^-1, "]")),
             col="cyan", lwd = 3)
        points(out350DF$nc, out350DF$NPP_VL, type="l", col="tomato", lwd = 3)
        points(equil350DF$nc_VL, equil350DF$NPP_VL, type="p", pch = 19, col = "blue")
        points(out350DF$nc, out350DF$NPP_350_L, type='l',col="violet", lwd = 3)
        points(out700DF$nc, out700DF$NPP_700, col="green", type="l", lwd = 3)
        points(equil350DF$nc_VL, inst700$equilNPP, type="p", col = "darkgreen", pch=19)
        points(equil700DF$nc_VL, equil700DF$NPP_VL, type="p", col="orange", pch = 19)
        points(equil700DF$nc_L, equil700DF$NPP_L,type="p", col="red", pch = 19)
        points(nfseq, NCMEDIUM$NPP, type="l", col="darkred", lwd = 3)
        points(Medium_equil_700$equilnf, Medium_equil_700$equilNPP, type="p", col="purple", pch = 19)
        
        
        legend("topright", c("P350", "P700", "VL", "L", "M",
                             "A", "B", "C", "D", "E"),
               col=c("cyan","green", "tomato", "violet","darkred","blue", "darkgreen","purple","red", "orange"), 
               lwd=c(2,2,2,2,2,NA,NA,NA,NA,NA), pch=c(NA,NA,NA,NA,NA,19,19,19,19,19), cex = 0.8, 
               bg = adjustcolor("grey", 0.8))
        
        
        dev.off()
        
    } else if (f.flag == 2) {
        return(cDF)
    } else if (f.flag == 3) {
        return(eDF)
    }
    
}
