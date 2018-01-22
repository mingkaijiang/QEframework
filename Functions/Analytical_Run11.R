
#### Analytical script to match GDAY Run 11 settings
####
#### Same as Run 10, except 
#### Considering medium term slow SOM pool without priming
#### i.e. only turning exudation on
################################################################################


#### Functions
Perform_Analytical_Run11 <- function(f.flag = 1, cDF, eDF) {
    #### Function to perform analytical run 9 simulations
    #### eDF: stores equilibrium points
    #### cDF: stores constraint points (curves)
    #### f.flag: = 1 simply plot analytical solution file
    #### f.flag: = 2 return cDF
    #### f.flag: = 3 return eDF

    ######### Main program
    source("Parameters/Analytical_Run10_Parameters.R")
    
    # create nc and pc for shoot to initiate
    nfseq <- round(seq(0.001, 0.1, by = 0.001),5)
    a_vec <- as.data.frame(allocn_exudation(nfseq))
    
    # plot photosynthetic constraints
    PC350 <- photo_constraint_full_cn(nfseq,a_vec,CO2=CO2_1)
    PC700 <- photo_constraint_full_cn(nfseq,a_vec,CO2=CO2_2)
    
    # calculate very long term NC and PC constraint on NPP, respectively
    NCVLONG <- NConsVLong(df=nfseq,a=a_vec)
    
    # solve very-long nutrient cycling constraint
    VLongN <- solveVLong_full_cn_medium(CO2_1)
    
    # Get Cpassive from very-long nutrient cycling solution
    aequiln <- allocn_exudation(VLongN$equilnf)
    pass <- slow_pool(df=VLongN$equilnf, a=aequiln)
    omegap <- aequiln$af*pass$omegafp + aequiln$ar*pass$omegarp
    CpassVLong <- omegap*VLongN$equilNPP/pass$decomp_p/(1-pass$qpq)*1000.0
    
    # Get Cslow from long nutrient cycling solution
    omegas <- aequiln$af*pass$omegafs + aequiln$ar*pass$omegars
    CslowLong <- omegas*VLongN$equilNPP/pass$decomp_s/(1-pass$qsq)*1000.0
    
    # Calculate nutrient release from woody pool
    NrelwoodVLong <- aequiln$aw*aequiln$nw*VLongN$equilNPP*1000.0
    
    # Calculate long term nutrient constraint
    NCHUGH <- NConsLong(nfseq, a_vec, CpassVLong,
                        NinL = Nin)
    
    # Calculate medium term nutrient constraint
    NCMEDIUM <- NConsMedium(df=nfseq, 
                            a=a_vec, 
                            Cpass=CpassVLong, 
                            Cslow=CslowLong, 
                            NinL = Nin+NrelwoodVLong)
    
    
    # Solve long-term equilibrium
    equil_long_350 <- solveLong_full_cn_medium(CO2=CO2_1, Cpass=CpassVLong, NinL = Nin)
    equil_long_700 <- solveLong_full_cn_medium(CO2=CO2_2, Cpass=CpassVLong, NinL = Nin)
    
    # Solve medium equilibrium
    equil_medium_350 <- solveMedium(CO2=CO2_1, Cpass=CpassVLong, Cslow=CslowLong, NinL = Nin+NrelwoodVLong)
    equil_medium_700 <- solveMedium(CO2=CO2_2, Cpass=CpassVLong, Cslow=CslowLong, NinL = Nin+NrelwoodVLong)
    
    # get the point instantaneous NPP response to doubling of CO2
    df700 <- as.data.frame(cbind(round(nfseq,3), PC700))
    inst700 <- inst_NPP(VLongN$equilnf, df700)
    
    ## locate the intersect between VL nutrient constraint and CO2 = 700
    VLong700 <- solveVLong_full_cn_medium(CO2_2)
    
    if (f.flag == 1) {
        
        ############################# Plotting
        ### plot 2-d plots of nf vs. npp and nf vs. pf
        tiff("Plots/Analytical_Run11_2d.tiff",
             width = 8, height = 7, units = "in", res = 300)
        par(mar=c(5.1,5.1,2.1,2.1))
        
        # Photosynthetic constraint CO2 = 350 ppm
        plot(nfseq,PC350,axes=T,
             type='l',xlim=c(0,0.05),ylim=c(0,3), 
             ylab = expression(paste("Production [kg C ", m^-2, " ", yr^-1, "]")),
             xlab = "Shoot N:C ratio", lwd = 2.5, col="cyan", cex.lab = 1.5)
        
        # Photosynthetic constraint CO2 = 700 ppm
        points(nfseq,PC700,type='l',col="green", lwd = 2.5)
        
        # VL nutrient constraint curve
        points(nfseq,NCVLONG$NPP_N,type='l',col="tomato", lwd = 2.5)
        
        # L nutrient constraint curve
        points(nfseq,NCHUGH$NPP,type='l',col="violet", lwd = 2.5)
        
        # M nutrient constraint curve
        points(nfseq,NCMEDIUM$NPP,type='l',col="darkred", lwd = 2.5)
        
        # VL intersect with CO2 = 350 ppm
        points(VLongN$equilnf,VLongN$equilNPP, pch = 19, cex = 2.0, col = "blue")
        
        # L intersect with CO2 = 700 ppm
        with(equil_long_700,points(equilnf,equilNPP,pch=19, cex = 2.0, col = "red"))
        
        # M intersect with CO2 = 700 ppm
        with(equil_medium_700,points(equilnf,equilNPP,pch=19, cex = 2.0, col = "purple"))
        
        # instantaneous NPP response to doubling CO2
        points(VLongN$equilnf, inst700$equilNPP, cex = 2.0, col = "darkgreen", pch=19)
        
        # VL intersect with CO2 = 700 ppm
        points(VLong700$equilnf, VLong700$equilNPP, cex = 2.0, col = "orange", pch = 19)
        
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