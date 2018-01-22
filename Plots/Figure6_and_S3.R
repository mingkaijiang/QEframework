#### To generate figure s3
#### explicit mineral N uptake as a coefficient
#### compare different uptake rates

#### programs

#Nutrient_uptake_comparison <- function() {
    
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
    
    ### new nuptakerate = 0.5
    nuptakerate <- 0.5  

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
    NCMEDIUM_2 <- NConsMedium_expl_min(nfseq, a_nf,CpassVLong, CslowLong,
                                     NinL = Nin+NrelwoodVLong)
    
    out350DF_2 <- data.frame(nfseq, NC350, NCVLONG, NCHUGH)
    colnames(out350DF_2) <- c("nc", "NPP_350", "NPP_VL",
                              "nleach_VL", "NPP_350_L", "nwood_L", "nburial_L",
                              "nleach_L", "aw")
    equil350DF_2 <- data.frame(VLongN, LongN)
    colnames(equil350DF_2) <- c("nc_VL", "NPP_VL", 
                                "nc_L", "NPP_L")
    
    ##### CO2 = 700
    # calculate NC vs. NPP at CO2 = 350 respectively
    NC700 <- photo_constraint_full_cn(nfseq, a_nf, CO2_2)
    
    # calculate very long term NC and PC constraint on NPP, respectively
    NCVLONG <- NConsVLong_expl_min(df=nfseq,a=a_nf)
    
    # solve very-long nutrient cycling constraint
    VLongN <- solveVLong_expl_min(CO2_2)
    
    out700DF_2 <- data.frame(nfseq, NC700, NCVLONG, NCHUGH)
    colnames(out700DF_2) <- c("nc", "NPP_700", "NPP_VL",
                              "nleach_VL", "NPP_700_L", "nwood_L", "nburial_L",
                              "nleach_L", "aw")
    
    # Find equilibrate intersection and plot
    LongN <- solveLong_expl_min(CO2_2, Cpass=CpassVLong, NinL=Nin)#+NrelwoodVLong)
    
    equil700DF_2 <- data.frame(VLongN, LongN)
    colnames(equil700DF_2) <- c("nc_VL", "NPP_VL", 
                                "nc_L", "NPP_L")
    
    # Find medium term equilibrium point
    Medium_equil_350_2 <- solveMedium_expl_min(CO2_1, Cpass = CpassVLong, Cslow = CslowLong, 
                                               NinL=Nin+NrelwoodVLong)
    Medium_equil_700_2 <- solveMedium_expl_min(CO2_2, Cpass = CpassVLong, Cslow = CslowLong, 
                                               NinL=Nin+NrelwoodVLong)
    
    # get the point instantaneous NPP response to doubling of CO2
    df700 <- as.data.frame(cbind(round(nfseq,3), NC700))
    inst700_2 <- inst_NPP(equil350DF_2$nc_VL, df700)
    
    #### nuptakerate = 1.5
    nuptakerate <- 1.5
    
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
    NCMEDIUM_3 <- NConsMedium_expl_min(nfseq, a_nf,CpassVLong, CslowLong,
                                     NinL = Nin+NrelwoodVLong)
    
    out350DF_3 <- data.frame(nfseq, NC350, NCVLONG, NCHUGH)
    colnames(out350DF_3) <- c("nc", "NPP_350", "NPP_VL",
                              "nleach_VL", "NPP_350_L", "nwood_L", "nburial_L",
                              "nleach_L", "aw")
    equil350DF_3 <- data.frame(VLongN, LongN)
    colnames(equil350DF_3) <- c("nc_VL", "NPP_VL", 
                                "nc_L", "NPP_L")
    
    ##### CO2 = 700
    # calculate NC vs. NPP at CO2 = 350 respectively
    NC700 <- photo_constraint_full_cn(nfseq, a_nf, CO2_2)
    
    # calculate very long term NC and PC constraint on NPP, respectively
    NCVLONG <- NConsVLong_expl_min(df=nfseq,a=a_nf)
    
    # solve very-long nutrient cycling constraint
    VLongN <- solveVLong_expl_min(CO2_2)
    
    out700DF_3 <- data.frame(nfseq, NC700, NCVLONG, NCHUGH)
    colnames(out700DF_3) <- c("nc", "NPP_700", "NPP_VL",
                              "nleach_VL", "NPP_700_L", "nwood_L", "nburial_L",
                              "nleach_L", "aw")
    
    # Find equilibrate intersection and plot
    LongN <- solveLong_expl_min(CO2_2, Cpass=CpassVLong, NinL=Nin)#+NrelwoodVLong)
    
    equil700DF_3 <- data.frame(VLongN, LongN)
    colnames(equil700DF_3) <- c("nc_VL", "NPP_VL", 
                                "nc_L", "NPP_L")
    
    # Find medium term equilibrium point
    Medium_equil_350_3 <- solveMedium_expl_min(CO2_1, Cpass = CpassVLong, Cslow = CslowLong, 
                                               NinL=Nin+NrelwoodVLong)
    Medium_equil_700_3 <- solveMedium_expl_min(CO2_2, Cpass = CpassVLong, Cslow = CslowLong, 
                                               NinL=Nin+NrelwoodVLong)
    
    # get the point instantaneous NPP response to doubling of CO2
    df700 <- as.data.frame(cbind(round(nfseq,3), NC700))
    inst700_3 <- inst_NPP(equil350DF_3$nc_VL, df700)
    
    
    #### Plotting
    ### plot 2-d plots of nf vs. npp 
    tiff("Plots/FigureS3.tiff",
         width = 12, height = 5, units = "in", res = 300)
    par(mfrow=c(1,3), mar=c(5.1,6.1,2.1,2.1))
    
    
    # shoot nc vs. NPP for nuptakerate = 0.5
    plot(out350DF_2$nc, out350DF_2$NPP_350, xlim=c(0.0, 0.05),
         ylim=c(0, 3), 
         type = "l", xlab = "Shoot N:C ratio", 
         ylab = expression(paste("Production [kg C ", m^-2, " ", yr^-1, "]")),
         col="cyan", lwd = 3)
    points(out350DF_2$nc, out350DF_2$NPP_VL, type="l", col="tomato", lwd = 3)
    points(equil350DF_2$nc_VL, equil350DF_2$NPP_VL, type="p", pch = 19, col = "blue", cex=1.5)
    points(out350DF_2$nc, out350DF_2$NPP_350_L, type='l',col="violet", lwd = 3)
    points(out700DF_2$nc, out700DF_2$NPP_700, col="green", type="l", lwd = 3)
    points(equil350DF_2$nc_VL, inst700_2$equilNPP, type="p", col = "darkgreen", pch=19, cex=1.5)
    points(equil700DF_2$nc_VL, equil700DF_2$NPP_VL, type="p", col="orange", pch = 19, cex = 1.5)
    points(equil700DF_2$nc_L, equil700DF_2$NPP_L,type="p", col="red", pch = 19, cex = 1.5)
    points(nfseq, NCMEDIUM_2$NPP, type="l", col="darkred", lwd = 3)
    points(Medium_equil_700_2$equilnf, Medium_equil_700_2$equilNPP, type="p", col="purple", pch = 19, cex = 1.5)
    
    # shoot nc vs. NPP for nuptakerate = 1.0
    plot(out350DF_1$nc, out350DF_1$NPP_350, xlim=c(0.0, 0.05),
         ylim=c(0, 3), 
         type = "l", xlab = "Shoot N:C ratio", 
         ylab = expression(paste("Production [kg C ", m^-2, " ", yr^-1, "]")),
         col="cyan", lwd = 3)
    points(out350DF_1$nc, out350DF_1$NPP_VL, type="l", col="tomato", lwd = 3)
    points(equil350DF_1$nc_VL, equil350DF_1$NPP_VL, type="p", pch = 19, col = "blue", cex = 1.5)
    points(out350DF_1$nc, out350DF_1$NPP_350_L, type='l',col="violet", lwd = 3)
    points(out700DF_1$nc, out700DF_1$NPP_700, col="green", type="l", lwd = 3)
    points(equil350DF_1$nc_VL, inst700_1$equilNPP, type="p", col = "darkgreen", pch=19, cex = 1.5)
    points(equil700DF_1$nc_VL, equil700DF_1$NPP_VL, type="p", col="orange", pch = 19, cex = 1.5)
    points(equil700DF_1$nc_L, equil700DF_1$NPP_L,type="p", col="red", pch = 19, cex = 1.5)
    points(nfseq, NCMEDIUM_1$NPP, type="l", col="darkred", lwd = 3)
    points(Medium_equil_700_1$equilnf, Medium_equil_700_1$equilNPP, type="p", col="purple", pch = 19, cex = 1.5)
    
    # shoot nc vs. NPP for nuptakerate = 1.5
    plot(out350DF_3$nc, out350DF_3$NPP_350, xlim=c(0.0, 0.05),
         ylim=c(0, 3), 
         type = "l", xlab = "Shoot N:C ratio", 
         ylab = expression(paste("Production [kg C ", m^-2, " ", yr^-1, "]")),
         col="cyan", lwd = 3)
    points(out350DF_3$nc, out350DF_3$NPP_VL, type="l", col="tomato", lwd = 3)
    points(equil350DF_3$nc_VL, equil350DF_3$NPP_VL, type="p", pch = 19, col = "blue", cex = 1.5)
    points(out350DF_3$nc, out350DF_3$NPP_350_L, type='l',col="violet", lwd = 3)
    points(out700DF_3$nc, out700DF_3$NPP_700, col="green", type="l", lwd = 3)
    points(equil350DF_3$nc_VL, inst700_3$equilNPP, type="p", col = "darkgreen", pch=19, cex = 1.5)
    points(equil700DF_3$nc_VL, equil700DF_3$NPP_VL, type="p", col="orange", pch = 19, cex = 1.5)
    points(equil700DF_3$nc_L, equil700DF_3$NPP_L,type="p", col="red", pch = 19, cex = 1.5)
    points(nfseq, NCMEDIUM_3$NPP, type="l", col="darkred", lwd = 3)
    points(Medium_equil_700_3$equilnf, Medium_equil_700_3$equilNPP, type="p", col="purple", pch = 19, cex = 1.5)
    
    
    legend("topright", c("P350", "P700", "VL", "L", "M",
                         "A", "B", "C", "D", "E"),
           col=c("cyan","green", "tomato", "violet","darkred","blue", "darkgreen","purple","red", "orange"), 
           lwd=c(2,2,2,2,2,NA,NA,NA,NA,NA), pch=c(NA,NA,NA,NA,NA,19,19,19,19,19), cex = 0.8, 
           bg = adjustcolor("grey", 0.8))
    
    
    dev.off()
    
    #### Plot Figure 6
    plotDF <- matrix(ncol=3, nrow = 12)
    colnames(plotDF) <- c("Model", "Timescale","Production")
    plotDF <- as.data.frame(plotDF)
    plotDF$Model <- rep(c("0.5", "1.0", "1.5"), each = 4)
    plotDF$Timescale <- rep(c("I", "M", "L", "VL"), 3)
    
    plotDF[plotDF$Model == "0.5" & plotDF$Timescale == "I", "Production"] <- (inst700_2$equilNPP - equil350DF_2$NPP_VL) / equil350DF_2$NPP_VL * 100
    plotDF[plotDF$Model == "0.5" & plotDF$Timescale == "M", "Production"] <- (Medium_equil_700_2$equilNPP - equil350DF_2$NPP_VL) / equil350DF_2$NPP_VL * 100
    plotDF[plotDF$Model == "0.5" & plotDF$Timescale == "L", "Production"] <- (equil700DF_2$NPP_L - equil350DF_2$NPP_VL) / equil350DF_2$NPP_VL * 100
    plotDF[plotDF$Model == "0.5" & plotDF$Timescale == "VL", "Production"] <-(equil700DF_2$NPP_VL - equil350DF_2$NPP_VL) / equil350DF_2$NPP_VL * 100
        
    plotDF[plotDF$Model == "1.0" & plotDF$Timescale == "I", "Production"] <- (inst700_1$equilNPP - equil350DF_1$NPP_VL) / equil350DF_1$NPP_VL * 100
    plotDF[plotDF$Model == "1.0" & plotDF$Timescale == "M", "Production"] <- (Medium_equil_700_1$equilNPP - equil350DF_1$NPP_VL) / equil350DF_1$NPP_VL * 100
    plotDF[plotDF$Model == "1.0" & plotDF$Timescale == "L", "Production"] <- (equil700DF_1$NPP_L - equil350DF_1$NPP_VL) / equil350DF_1$NPP_VL * 100
    plotDF[plotDF$Model == "1.0" & plotDF$Timescale == "VL", "Production"] <- (equil700DF_1$NPP_VL - equil350DF_1$NPP_VL) / equil350DF_1$NPP_VL * 100
    
    plotDF[plotDF$Model == "1.5" & plotDF$Timescale == "I", "Production"] <- (inst700_3$equilNPP - equil350DF_3$NPP_VL) / equil350DF_3$NPP_VL * 100
    plotDF[plotDF$Model == "1.5" & plotDF$Timescale == "M", "Production"] <- (Medium_equil_700_3$equilNPP - equil350DF_3$NPP_VL) / equil350DF_3$NPP_VL * 100
    plotDF[plotDF$Model == "1.5" & plotDF$Timescale == "L", "Production"] <- (equil700DF_3$NPP_L - equil350DF_3$NPP_VL) / equil350DF_3$NPP_VL * 100 
    plotDF[plotDF$Model == "1.5" & plotDF$Timescale == "VL", "Production"] <-(equil700DF_3$NPP_VL - equil350DF_3$NPP_VL) / equil350DF_3$NPP_VL * 100
    
    plotDF$Timescale <- factor(plotDF$Timescale, levels=c("I", "M", "L", "VL"))
    
    require(ggplot2)
    
    ylabel <- bquote(.("Response to") ~ eCO[2] ~ ("%"))
    
    
    # making bar plots
    tiff("Plots/Figure6.tiff",
         width = 10, height = 5, units = "in", res = 300)

    p1 <- ggplot(plotDF, aes(x=Model, y=Production, fill=Timescale)) +   
        geom_bar(position='dodge', stat='identity') +
        labs(list(x = "N uptake rate", y = ylabel, fill = "Timescale")) +
        scale_y_continuous(limits=c(-10,30)) +
        scale_fill_manual(values=c("darkgreen", "purple", "red", "orange")) +
        theme(text = element_text(size=20),
              axis.text.x = element_text(),
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              axis.line = element_line(colour = "black")) 
    print(p1)
    
    dev.off()
    
#}




#### Script
#Nutrient_uptake_comparison()
