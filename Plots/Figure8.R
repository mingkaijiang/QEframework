#### To plot Figure 8
#### comparison of analytical run 10 and 11
#### Run 10 is exudation on, priming on
#### Run 11 is exudation on, priming off

#### Program

priming_effect_plot <- function() {

    #### The parameters are the same
    source("Parameters/Analytical_Run10_Parameters.R")
    
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
    equil_long_350_on <- solveLong_full_cn_medium(CO2=CO2_1, Cpass=CpassVLong, NinL = Nin)
    equil_long_700_on <- solveLong_full_cn_medium(CO2=CO2_2, Cpass=CpassVLong, NinL = Nin)
    
    # Solve medium equilibrium
    equil_medium_350_on <- solveMedium_priming(CO2=CO2_1, Cpass=CpassVLong, Cslow=CslowLong, NinL = Nin+NrelwoodVLong)
    equil_medium_700_on <- solveMedium_priming(CO2=CO2_2, Cpass=CpassVLong, Cslow=CslowLong, NinL = Nin+NrelwoodVLong)
    
    # get the point instantaneous NPP response to doubling of CO2
    df700 <- as.data.frame(cbind(round(nfseq_on,3), PC700_on))
    inst700_on <- inst_NPP(VLongN_on$equilnf, df700)
    
    ## locate the intersect between VL nutrient constraint and CO2 = 700
    VLong700_on <- solveVLong_exudation_medium(CO2_2)
    
    
    ######## Analytical run 11 - exudation on, priming off
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
    equil_long_350_off <- solveLong_full_cn_medium(CO2=CO2_1, Cpass=CpassVLong, NinL = Nin)#+NrelwoodVLong)
    equil_long_700_off <- solveLong_full_cn_medium(CO2=CO2_2, Cpass=CpassVLong, NinL = Nin)#+NrelwoodVLong)
    
    # Solve medium equilibrium
    equil_medium_350_off <- solveMedium(CO2=CO2_1, Cpass=CpassVLong, Cslow=CslowLong, NinL = Nin+NrelwoodVLong)
    equil_medium_700_off <- solveMedium(CO2=CO2_2, Cpass=CpassVLong, Cslow=CslowLong, NinL = Nin+NrelwoodVLong)
    
    # get the point instantaneous NPP response to doubling of CO2
    df700 <- as.data.frame(cbind(round(nfseq_off,3), PC700_off))
    inst700_off <- inst_NPP(VLongN_off$equilnf, df700)
    
    ## locate the intersect between VL nutrient constraint and CO2 = 700
    VLong700_off <- solveVLong_full_cn_medium(CO2_2)
    
    
    ########## Plotting
    tiff("Plots/Figure8.tiff",
         width = 12, height = 7, units = "in", res = 300)
    par(mfrow=c(1,2), mar=c(5.1,6.1,2.1,2.1))
    
    
    # priming off
    plot(nfseq_off,PC350_off,axes=T,
         type='l',xlim=c(0,0.02),ylim=c(1,2), 
         ylab = expression(paste("Production [kg C ", m^-2, " ", yr^-1, "]")),
         xlab = "Shoot N:C ratio", lwd = 2.5, col="cyan", cex.lab = 1.5)
    points(nfseq_off,PC700_off,type='l',col="green", lwd = 2.5)
    points(nfseq_off,NCVLONG_off$NPP_N,type='l',col="tomato", lwd = 2.5)
    points(nfseq_off,NCHUGH_off$NPP,type='l',col="violet", lwd = 2.5)
    points(nfseq_off,NCMEDIUM_off$NPP,type='l',col="darkred", lwd = 2.5)
    points(VLongN_off$equilnf,VLongN_off$equilNPP, pch = 19, cex = 2.0, col = "blue")
    with(equil_long_700_off,points(equilnf,equilNPP,pch=19, cex = 2.0, col = "red"))
    with(equil_medium_700_off,points(equilnf,equilNPP,pch=19, cex = 2.0, col = "purple"))
    points(VLongN_off$equilnf, inst700_off$equilNPP, cex = 2.0, col = "darkgreen", pch=19)
    points(VLong700_off$equilnf, VLong700_off$equilNPP, cex = 2.0, col = "orange", pch = 19)
    text(x=0.001, y=2.9, "(a)", cex = 2)
    
    
    # priming on
    plot(nfseq_on,PC350_on,axes=T,
         type='l',xlim=c(0,0.02),ylim=c(1,2), 
         ylab = expression(paste("Production [kg C ", m^-2, " ", yr^-1, "]")),
         xlab = "Shoot N:C ratio", lwd = 2.5, col="cyan", cex.lab = 1.5)
    points(nfseq_on,PC700_on,type='l',col="green", lwd = 2.5)
    points(nfseq_on,NCVLONG_on$NPP_N,type='l',col="tomato", lwd = 2.5)
    points(nfseq_on,NCHUGH_on$NPP,type='l',col="violet", lwd = 2.5)
    points(nfseq_on,NCMEDIUM_on$NPP,type='l',col="darkred", lwd = 2.5)
    points(VLongN_on$equilnf,VLongN_on$equilNPP, pch = 19, cex = 2.0, col = "blue")
    with(equil_long_700_on,points(equilnf,equilNPP,pch=19, cex = 2.0, col = "red"))
    with(equil_medium_700_on,points(equilnf,equilNPP,pch=19, cex = 2.0, col = "purple"))
    points(VLongN_on$equilnf, inst700_on$equilNPP, cex = 2.0, col = "darkgreen", pch=19)
    points(VLong700_on$equilnf, VLong700_on$equilNPP, cex = 2.0, col = "orange", pch = 19)
    text(x=0.001, y=2.9, "(b)", cex = 2)
    
    
    # legend
    legend("topright", c("P350", "P700", "VL", "L", "M",
                         "A", "B", "C", "D", "E"),
           col=c("cyan","green", "tomato", "violet","darkred","blue", "darkgreen","purple","red", "orange"), 
           lwd=c(2,2,2,2,2,NA,NA,NA,NA,NA), pch=c(NA,NA,NA,NA,NA,19,19,19,19,19), cex = 1.0, 
           bg = adjustcolor("grey", 0.8), ncol=2)
    
    dev.off()
    
}



#### Script
priming_effect_plot()

