


Static_Simulation_Analytical_Plot <- function() {
    
    # library
    require(plotly)
    
    ### Create df to store all the constraints
    nfseq <- round(seq(0.001, 0.1, b=0.001), 5)
    csDF <- data.frame(nfseq, NA, NA, NA, NA, NA)
    colnames(csDF) <- c("nf", "aCO2", "eCO2", "VL",
                        "L", "M")
    
    
    ######### Main program
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
    # PCMEDIUM_350 is implicit, but can also be calculated if needed
    
    Medium_equil_350 <- solveMedium_full_cnp(CO2=CO2_1, Cpass=CpassVLong, Cslow=CslowLong, NinL = Nin+NrelwoodVLong,
                                             PinL=Pin+PrelwoodVLong)
    
    csDF$aCO2 <- Photo350
    csDF$VL <- NCVLONG$NPP_N
    csDF$L <- NCLONG$NPP
    csDF$M <- NCMEDIUM$NPP
    
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
    
    csDF$eCO2 <- Photo700
    
    
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
    
    
    ### time series data
    mvDF <- read.csv("GDAY/analyses/Run1/annual_gday_result_transient_CO2_ELE.csv", header=T, sep=",")
    mvDF$nc <- mvDF$shootn/mvDF$shoot
    mvDF$npp_new <- mvDF$npp/10
    with(mvDF[1:100, ], plot(npp_new~year))
    
    ### Plotting
    tiff("Plots/static_analytical_comparison.tiff",
         width = 6.5, height = 6.5, units = "in", res = 200)
    par(mar=c(5.1,6.1,2.1,2.1))
    
    split.screen(c(2,2))
    
    screen(3)
    # plot the baseline constraint curves
    with(csDF, plot(aCO2~nf, type="l", 
                        xlim=c(0.01, 0.03),
                        ylim=c(0.8, 2.5), 
                        xlab = "Shoot N:C ratio", 
                        ylab = expression(paste("Production [kg C ", m^-2, " ", yr^-1, "]")),
                        col="cyan", lwd = 1.5, cex.lab=1.0))
    with(csDF, lines(eCO2~nf, col="green", type="l", lwd = 1.5))
    with(csDF, lines(VL~nf, type="l", col="tomato", lwd = 1.5))
    with(csDF, lines(L~nf, type="l", col="violet", lwd = 1.5))
    with(csDF, lines(M~nf, type="l", col="darkred", lwd = 1.5))
    points(inst700$nf, inst700$equilNPP, cex=2, pch=16, col="darkgreen")
    points(equil350DF$nc_VL, equil350DF$NPP_VL, cex=2, pch=16, col="blue")
    points(equil700DF$nc_L, equil700DF$NPP_L, cex=2, pch=16, col="red")
    points(equil700DF$nc_VL, equil700DF$NPP_VL, cex=2, pch=16, col="orange")
    points(Medium_equil_700$equilnf, Medium_equil_700$equilNPP, cex=2, pch=16, col="purple")
    
    
    # plot the time series nf
    screen(1)
    with(mvDF[10:5000,], plot(year, nc+0.003, pch=1, cex=0.1,
                         xlim=c(0, 600), ylim=c(0.01, 0.03),
                         xlab="Year", ylab="Shoot N:C ratio", cex.lab=1.0, col="black"))
    points(0, mvDF$nc[5], col="blue", cex=2.0, pch=16)
    with(mvDF[6, ], points(year, nc, col="darkgreen", cex=2.0, pch=16))
    with(mvDF[26, ], points(year, nc+0.003, col="purple", cex=2.0, pch=16))
    with(mvDF[206, ], points(year, nc+0.003, col="red", cex=2.0, pch=16))
    with(mvDF[625, ], points(year, nc+0.003, col="orange", cex=2.0, pch=16))
    
    
    # plot time series NPP
    screen(2)
    with(mvDF[7:5000,], plot(year, npp_new+0.05, pch=1, cex=0.1, type="p",
                              xlim=c(0, 600), ylim=c(1.4, 2),
                              xlab="Year", ylab=expression(paste("Production [kg C ", m^-2, " ", yr^-1, "]")),
                              cex.lab=1.0, col="black"))
    with(mvDF[5, ], points(year, npp_new+0.05, col="blue", cex=2.0, pch=16))
    with(mvDF[6, ], points(year, npp_new+0.15, col="darkgreen", cex=2.0, pch=16))
    with(mvDF[26, ], points(year, npp_new+0.05, col="purple", cex=2.0, pch=16))
    with(mvDF[206, ], points(year, npp_new+0.05, col="red", cex=2.0, pch=16))
    with(mvDF[625, ], points(year, npp_new+0.05, col="orange", cex=2.0, pch=16))
    points(mvDF$year[6], mvDF$npp_new[8]+0.15, pch=1, cex=0.1, col="black")
    points(mvDF$year[6], mvDF$npp_new[9]+0.15, pch=1, cex=0.1, col="black")
    points(mvDF$year[6], mvDF$npp_new[10]+0.15, pch=1, cex=0.1, col="black")
    points(mvDF$year[6], mvDF$npp_new[11]+0.15, pch=1, cex=0.1, col="black")
    points(mvDF$year[6], mvDF$npp_new[12]+0.15, pch=1, cex=0.1, col="black")
    
    screen(4)
    plot(1,xaxt='n',yaxt='n',bty='n',pch='',ylab='',xlab='')
    # add legends
    legend("bottomright", c("P350", "P700", "VL", "L", "M",
                            "A", "B", "C", "D", "E"),
           col=c("cyan","green", "tomato", "violet","darkred","blue", "darkgreen","purple","red", "orange"), 
           lwd=c(2,2,2,2,2,NA,NA,NA,NA,NA), pch=c(NA,NA,NA,NA,NA,19,19,19,19,19), cex = 1.0, 
           bg = adjustcolor("grey", 0.8), ncol=2)
    
    
    dev.off()
    
}
