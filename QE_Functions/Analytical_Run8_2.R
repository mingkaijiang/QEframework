
#### Analytical script to match GDAY Run 8 settings
####
#### Same as Run 8.1, except
#### 1. N uptake rates as a function of root biomass - O-CN approach: saturaing function of mineral N
#### 2. Fixed passive SOM stoichiometry
#### 3. Using N only model
################################################################################


#### Functions
Perform_Analytical_Run8_2 <- function() {
    #### Function to perform analytical run 8.2 simulations

    ######### Main program
    source("Parameters/Analytical_Run8_2_Parameters.R")
    
    # N:C ratios for x-axis
    nfseq <- seq(0.001,0.1,by=0.001)
    # need allocation fractions here
    a_vec <- allocn(nfseq)
    
    # plot photosynthetic constraints
    PC350 <- photo_constraint_full_cn(nfseq,a_vec,CO2=CO2_1)
    PC700 <- photo_constraint_full_cn(nfseq,a_vec,CO2=CO2_2)
    
    # calculate very long term NC and PC constraint on NPP, respectively
    VLongN <- NConsVLong_root_ocn(CO2_1)

    # Get Cpassive from very-long nutrient cycling solution
    aequiln <- allocn(VLongN$equilnf)
    pass <- slow_pool(df=VLongN$equilnf, a=aequiln)
    omegap <- aequiln$af*pass$omegafp + aequiln$ar*pass$omegarp
    omegas <- aequiln$af*pass$omegafs + aequiln$ar*pass$omegars
    CpassVLong <- omegap*VLongN$equilNPP/pass$decomp_p/(1-pass$qpq)*1000.0
    
    # Calculate nutrient release from recalcitrant pools
    NrelwoodVLong <- aequiln$aw*aequiln$nw*VLongN$equilNPP*1000.0
    
    # Calculate long term nutrieng constraint
    NCHUGH <- NConsLong_root_ocn(df=nfseq, a=a_vec,Cpass=CpassVLong,
                                  NinL = Nin)#+NrelwoodVLong)
    
    # Find equilibrate intersection and plot
    equil_long_350 <- solveLong_root_ocn(CO2_1, Cpass=CpassVLong, NinL= Nin)#+NrelwoodVLong)
    equil_long_700 <- solveLong_root_ocn(CO2_2, Cpass=CpassVLong, NinL= Nin)#+NrelwoodVLong)
    
    CslowLong <- omegas*equil_long_350$equilNPP/pass$decomp_s/(1-pass$qsq)*1000.0
    
    # plot medium nutrient cycling constraint
    NCMEDIUM <- NConsMedium_root_ocn(nfseq, a_vec, Cpass=CpassVLong, Cslow=CslowLong, NinL=Nin+NrelwoodVLong)
    
    # solve medium term equilibrium at CO2 = 700 ppm
    equil_medium_700 <- solveMedium_root_ocn(CO2_2,Cpass=CpassVLong,Cslow=CslowLong,Nin=Nin+NrelwoodVLong)
    

    # get the point instantaneous NPP response to doubling of CO2
    df700 <- as.data.frame(cbind(round(nfseq,3), PC700))
    inst700 <- inst_NPP(VLongN$equilnf, df700)
    
    ## locate the intersect between VL nutrient constraint and CO2 = 700
    VLong700 <-  NConsVLong_root_ocn(CO2_2)
    
    
    ######### Plotting
    
    tiff("Plots/Analytical_Run8_2.tiff",
         width = 8, height = 7, units = "in", res = 300)
    par(mar=c(5.1,5.1,2.1,2.1))

    
    # shoot nc vs. NPP
    # Photosynthetic constraint CO2 = 350 ppm
    plot(nfseq,PC350,axes=T,
         type='l',xlim=c(0,0.05),ylim=c(0,3), 
         ylab = expression(paste("Production [kg C ", m^-2, " ", yr^-1, "]")),
         xlab = "Shoot N:C ratio", lwd = 2.5, col="cyan", cex.lab = 1.5)

    # Photosynthetic constraint CO2 = 700 ppm
    points(nfseq,PC700,type='l',col="green", lwd = 2.5)
    
    # L nutrient constraint curve
    points(nfseq,NCHUGH$NPP,type='l',col="violet", lwd = 2.5)
    
    # VL intersect with CO2 = 350 ppm
    points(VLongN$equilnf,VLongN$equilNPP, pch = 19, cex = 2.0, col = "blue")
    
    # L intersect with CO2 = 700 ppm
    with(equil_long_700,points(equilnf,equilNPP,pch=19, cex = 2.0, col = "red"))
    
    
    # VL intersect with CO2 = 700 ppm
    points(VLong700$equilnf, VLong700$equilNPP, cex = 2.0, col = "orange", pch = 19)
    
    # M nutrient curve
    points(nfseq, NCMEDIUM$NPP, type="l", col="darkred", lwd = 2.5)
    
    # M intersect with CO2 = 700 ppm
    points(equil_medium_700$equilnf, equil_medium_700$equilNPP, cex = 2.0, col = "purple", pch = 19)
    
    # instantaneous NPP response to doubling CO2
    points(VLongN$equilnf, inst700$equilNPP, cex = 1, col = "darkgreen", pch=19)
    
    legend("topright", c("P350", "P700", "VL", "L", "M",
                         "A", "B", "C", "D", "E"),
           col=c("cyan","green", "tomato", "violet","darkred","blue", "darkgreen","purple","red", "orange"), 
           lwd=c(2,2,2,2,2,NA,NA,NA,NA,NA), pch=c(NA,NA,NA,NA,NA,19,19,19,19,19), cex = 0.8, 
           bg = adjustcolor("grey", 0.8))
    
    
    dev.off()    

    
}
