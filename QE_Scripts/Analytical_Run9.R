
#### Analytical script to match GDAY Run 9 settings
####
#### Same as Run 7, except 
#### 1. turn root exudation on
#### 2. turn priming effect on (effect on passive SOM pool)
#### 3. N only model (same as Run 2, but with explicit mineral N pool)
#### 4. Fix wood stoichiometry (and therefore adjust wood nc parameters)
####
################################################################################


#### Functions
Perform_Analytical_Run9 <- function(f.flag = 1, cDF, eDF) {
    #### Function to perform analytical run 9 simulations
    #### eDF: stores equilibrium points
    #### cDF: stores constraint points (curves)
    #### f.flag: = 1 simply plot analytical solution file
    #### f.flag: = 2 return cDF
    #### f.flag: = 3 return eDF
    
    ######### Main program
    #### Function to perform analytical run 1 simulations
    #### eDF: stores equilibrium points
    #### cDF: stores constraint points (curves)
    #### f.flag: = 1 simply plot analytical solution file
    #### f.flag: = 2 return cDF
    #### f.flag: = 3 return eDF
    
    ######### Main program
    source("Parameters/Analytical_Run9_Parameters.R")
    
    # create nc and pc for shoot to initiate
    nfseq <- round(seq(0.01, 0.1, by = 0.001),5)
    a_vec <- as.data.frame(allocn_exudation(nfseq))
    
    # plot photosynthetic constraints
    PC350 <- photo_constraint_full_cn(nfseq,a_vec,CO2=CO2_1)
    PC700 <- photo_constraint_full_cn(nfseq,a_vec,CO2=CO2_2)
    
    # calculate very long term NC and PC constraint on NPP, respectively
    NCVLONG <- NConsVLong(df=nfseq,a=a_vec)

    # solve very-long nutrient cycling constraint
    VLongN <- solveVLong_exudation(CO2_1)

    # Get Cpassive from very-long nutrient cycling solution
    aequiln <- allocn_exudation(VLongN$equilnf)
    pass <- passive_exudation(df=VLongN$equilnf, a=aequiln, npp=VLongN$equilNPP_N)
    omega <- aequiln$af*pass$omegaf + aequiln$ar*pass$omegar
    CpassVLong <- omega*VLongN$equilNPP/pass$decomp/(1-pass$qq)*1000.0
    
    # Calculate nutrient release from recalcitrant pools
    NrelwoodVLong <- aequiln$aw*aequiln$nw*VLongN$equilNPP_N*1000.0
    
    # Calculate long term nutrieng constraint
    NCHUGH <- NConsLong_exudation(nfseq, a_vec, CpassVLong,
                                  NinL = Nin)#+NrelwoodVLong)
    
    # compute different passive pool based on original function
    pass_orig <- passive(df=VLongN$equilnf, a=aequiln)
    omega_orig <- aequiln$af*pass_orig$omegaf + aequiln$ar*pass_orig$omegar
    CpassVLong_orig <- omega_orig*VLongN$equilNPP/pass_orig$decomp/(1-pass_orig$qq)*1000.0
    NrelwoodVLong_orig <- aequiln$aw*aequiln$nw*VLongN$equilNPP_N*1000.0
    NCHUGH_orig <- NConsLong(nfseq, a_vec, CpassVLong_orig,
                                       NinL = Nin)#+NrelwoodVLong_orig)
    
    # Solve longterm equilibrium
    equil_long_350 <- solveLong_exudation(CO2=CO2_1, Cpass=CpassVLong, NinL = Nin)#+NrelwoodVLong)
    equil_long_700 <- solveLong_exudation(CO2=CO2_2, Cpass=CpassVLong, NinL = Nin)#+NrelwoodVLong)
    
    # get the point instantaneous NPP response to doubling of CO2
    df700 <- as.data.frame(cbind(round(nfseq,3), PC700))
    inst700 <- inst_NPP(VLongN$equilnf, df700)
    
    ## locate the intersect between VL nutrient constraint and CO2 = 700
    VLong700 <- solveVLong_exudation(CO2_2)
    

    if (f.flag == 1) {
        
        ######### Plotting
        ### plot 2-d plots of nf vs. npp and nf vs. pf
        tiff("Plots/Analytical_Run9_2d.tiff",
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
        
        points(nfseq,NCHUGH_orig$NPP,type='l',col="grey", lwd = 2.5)
        
        # VL intersect with CO2 = 350 ppm
        points(VLongN$equilnf,VLongN$equilNPP, pch = 19, cex = 2.0, col = "blue")
        
        # L intersect with CO2 = 700 ppm
        with(equil_long_700,points(equilnf,equilNPP,pch=19, cex = 2.0, col = "red"))
        
        # instantaneous NPP response to doubling CO2
        points(VLongN$equilnf, inst700$equilNPP, cex = 2.0, col = "darkgreen", pch=19)
        
        # VL intersect with CO2 = 700 ppm
        points(VLong700$equilnf, VLong700$equilNPP, cex = 2.0, col = "orange", pch = 19)
        
        #legend("topright", c(expression(paste("Photo constraint at ", CO[2]," = 350 ppm")), 
        #                     expression(paste("Photo constraint at ", CO[2]," = 700 ppm")), 
        #                     "VL nutrient constraint", "L nutrient constraint priming on",
        #                     "L nutrient constraint priming off",
        #                     "A", "B", "C", "D"),
        #       col=c("cyan","green", "tomato", "violet","grey","blue", "darkgreen", "red", "orange"), 
        #       lwd=c(2,2,2,2,2,NA,NA,NA,NA), pch=c(NA,NA,NA,NA,NA,19,19,19,19), cex = 1.0, 
        #       bg = adjustcolor("grey", 0.8))
        
        legend("topright", c("P350", "P700", "VL", "L priming on", "L priming off",
                             "A", "B", "C", "D"),
               col=c("cyan","green", "tomato", "violet", "grey", "blue", "darkgreen", "red", "orange"), 
               lwd=c(2,2,2,2,2,NA,NA,NA,NA), pch=c(NA,NA,NA,NA,NA,19,19,19,19), cex = 0.8, 
               bg = adjustcolor("grey", 0.8))
        
        
        dev.off()
        
    } else if (f.flag == 2) {
        return(cDF)
    } else if (f.flag == 3) {
        return(eDF)
    }
    
}