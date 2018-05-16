
#### Analytical script Run 8.5
####
#### Assumptions:
#### 1. turn root exudation on
#### 2. turn priming effect on (effect on passive SOM pool)
#### 3. N only model with explicit mineral N considered
#### 4. Variable wood stoichiometry 
#### 5. fcue = 0.5
####
################################################################################
#### Functions
Perform_Analytical_Run85 <- function(f.flag = 1) {
    #### Function to perform analytical run 85 simulations
    #### eDF: stores equilibrium points
    #### cDF: stores constraint points (curves)
    #### f.flag: = 1 simply plot analytical solution and create individual pdf file
    #### f.flag: = 2 return a list consisting of two dataframes

    ######### Main program
    source("Parameters/Analytical_Run85_Parameters.R")
    
    ### create a range of nc for shoot to initiate
    nfseq <- round(seq(0.001, 0.1, by = 0.001),5)
    
    ### create nc ratio for wood, root, and allocation coefficients
    a_nf <- as.data.frame(alloc_prim(nfseq))
    
    ### calculate photosynthetic constraint at CO2 = 350
    P350 <- photo_constraint_full(nf=nfseq, nfdf=a_nf, CO2=CO2_1)

    ### calculate very long term NC constraint on NPP, respectively
    VL <- VL_constraint(nf=nfseq, nfdf=a_nf)
    
    ### finding the equilibrium point between photosynthesis and very long term nutrient constraints
    VL_eq <- solve_VL_full_prim(CO2=CO2_1)

    ### calculate nw and nr for VL equilibrated nf value
    a_eq <- alloc_prim(VL_eq$nf)
    
    ### calculate soil parameters, e.g. reburial coef.
    s_coef <- soil_coef(df=VL_eq$nf, a=a_eq)
    
    ### note: a fraction of root is used for exudation
    ###       so they should be excluded here. 
    omega_ap <- a_eq$af*s_coef$omega_af_pass + (a_eq$ar-a_eq$ar*a_eq$ariz)*s_coef$omega_ar_pass
    omega_as <- a_eq$af*s_coef$omega_af_slow + (a_eq$ar-a_eq$ar*a_eq$ariz)*s_coef$omega_ar_slow
    
    ### Get C from very-long term nutrient cycling solution
    ### return in g C m-2 
    C_pass_VL <- omega_ap*VL_eq$NPP/s_coef$decomp_pass/(1-s_coef$qq_pass)*1000.0
    
    ### Calculate nutrient release from woody pool
    ### return in g N m-2 yr-1
    N_wood_L <- a_eq$aw*a_eq$nw*VL_eq$NPP*1000.0

    ### Calculate long term nutrient constraint
    L <- L_constraint_prim(df=nfseq, a=a_nf, 
                            C_pass=C_pass_VL,
                            Nin_L = Nin)
    
    ### Find long term equilibrium point
    L_eq <- solve_L_full_prim(CO2=CO2_1, C_pass=C_pass_VL, Nin_L = Nin)
    
    ### Get Cslow from long nutrient cycling solution
    ### return in g C m-2
    ### does not consider priming
    C_slow_L_prim_off <- omega_as*L_eq$NPP/s_coef$decomp_slow/(1-s_coef$qq_slow)*1000.0
    
    # calculate N and C gaps based on priming effect
    c_into_active <- VL_eq$NPP * a_eq$ar * a_eq$ariz * rhizo_cue * 1000.0
    n_into_active <- c_into_active * a_eq$nr
    n_active_gap <- c_into_active * nc_active - n_into_active
    
    # adjust decomposition of slow pool to close the N gap
    decomp_slow_prim_on <- s_coef$decomp_slow * (1 + km) * pmax(c_into_active/(c_into_active + km), 0.3)
    
    # Calculate C slow based on exudation and new decomposition values
    C_slow_L_prim_on <- omega_as*VL_eq$NPP/decomp_slow_prim_on/(1-s_coef$qq_slow)*1000.0
    
    # 2nd method to calculate new C slow pool
    C_slow_L_prim_on_2 <- C_slow_L_prim_off - (n_active_gap / ncs) 
    
    ### Calculate medium term nutrient constraint
    M <- M_constraint_prim(df=nfseq,a=a_nf, 
                           C_pass=C_pass_VL, 
                           C_slow=C_slow_L_prim_on, 
                           Nin_L = Nin+N_wood_L)
    
    ### calculate M equilibrium point
    M_eq <- solve_M_full_prim(CO2=CO2_1, 
                              C_pass=C_pass_VL, 
                              C_slow=C_slow_L_prim_on, 
                              Nin_L = Nin+N_wood_L)
    

    out350DF <- data.frame(CO2_1, nfseq, P350, VL$NPP, 
                           L$NPP, M$NPP)
    colnames(out350DF) <- c("CO2", "nc", "NPP_photo", "NPP_VL",
                            "NPP_L", "NPP_M")
    equil350DF <- data.frame(CO2_1, VL_eq, L_eq, M_eq)
    colnames(equil350DF) <- c("CO2", "nc_VL", "NPP_VL", 
                              "nc_L", "NPP_L", "nc_M", "NPP_M")
    
    ##### CO2 = 700
    ### photo constraint
    P700 <- photo_constraint_full(nf=nfseq, nfdf=a_nf, CO2=CO2_2)
    
    ### VL equilibrated point with eCO2
    VL_eq <- solve_VL_full_prim(CO2=CO2_2)
    
    ### Find long term equilibrium point
    L_eq <- solve_L_full_prim(CO2=CO2_2, C_pass=C_pass_VL, Nin_L = Nin)
    
    ### Find medium term equilibrium point
    M_eq <- solve_M_full_prim(CO2=CO2_2, 
                         C_pass=C_pass_VL, 
                         C_slow=C_slow_L_prim_on, 
                         Nin_L = Nin+N_wood_L)
    
    out700DF <- data.frame(CO2_2, nfseq, P700, 
                           VL$NPP, L$NPP, M$NPP)
    colnames(out700DF) <- c("CO2", "nc", "NPP_photo", "NPP_VL",
                            "NPP_L", "NPP_M")
    
    equil700DF <- data.frame(CO2_2, VL_eq, L_eq, M_eq)
    colnames(equil700DF) <- c("CO2", "nc_VL", "NPP_VL", 
                              "nc_L", "NPP_L", "nc_M", "NPP_M")
 
    ### get the point instantaneous NPP response to doubling of CO2
    df700 <- as.data.frame(cbind(round(nfseq,3), P700))
    inst700 <- inst_NPP(equil350DF$nc_VL, df700)
    
    equil350DF$NPP_I <- inst700$equilNPP
    equil700DF$NPP_I <- inst700$equilNPP
    
    if (f.flag == 1) {
        
        ### plot 2-d plots of nf vs. npp and nf vs. pf
        tiff("Plots/Analytical_Run85_2d.tiff",
             width = 5, height = 5, units = "in", res = 300)
        par(mar=c(5.1,6.1,2.1,2.1))
        
        ### shoot nc vs. NPP
        plot(out350DF$nc, out350DF$NPP_photo, xlim=c(0.001, 0.03),
              ylim=c(0, 2.0), 
             type = "l", xlab = "Leaf N:C ratio", 
             ylab = expression(paste("NPP [kg C ", m^-2, " ", yr^-1, "]")),
             col="cyan", lwd = 3, cex.lab=1.5)

        abline(h = seq(0.5, 3.0, 0.5), v = seq(0.01, 0.05, 0.01), col="lightgray", lty = 3)
        
        points(out350DF$nc, out350DF$NPP_VL, type="l", col="tomato", lwd = 3)
        
        points(equil350DF$nc_VL, equil350DF$NPP_VL, type="p", pch = 19, col = "blue", cex = 2)
        
        points(out350DF$nc, out350DF$NPP_L, type='l',col="violet", lwd = 3)
        
        points(out350DF$nc, out350DF$NPP_M, type="l", col="darkred", lwd = 3)
        
        points(out700DF$nc, out700DF$NPP_photo, col="green", type="l", lwd = 3)
        
        points(equil350DF$nc_VL, inst700$equilNPP, type="p", col = "darkgreen", pch=19, cex = 2)
        
        points(equil700DF$nc_VL, equil700DF$NPP_VL, type="p", col="orange", pch = 19, cex = 2)
        
        points(equil700DF$nc_L, equil700DF$NPP_L,type="p", col="red", pch = 19, cex = 2)
        
        points(equil700DF$nc_M, equil700DF$NPP_M, type="p", col="purple", pch = 19, cex = 2)
        
        legend("bottomright", c("P350", "P700", "VL", "L", "M",
                            "A", "B", "C", "D", "E"),
               col=c("cyan","green", "tomato", "violet","darkred","blue", "darkgreen","purple","red", "orange"), 
               lwd=c(2,2,2,2,2,NA,NA,NA,NA,NA), pch=c(NA,NA,NA,NA,NA,19,19,19,19,19), cex = 1.0, 
               bg = adjustcolor("grey", 0.8), ncol=2)
        
        dev.off()
        
    } else if (f.flag == 2) {
        
        my.list <- list(cDF = data.frame(rbind(out350DF, out700DF)), 
                        eDF = data.frame(rbind(equil350DF, equil700DF)))
        
        return(my.list)
    } 
}

