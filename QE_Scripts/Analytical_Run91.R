
#### Analytical script Run 9.1
####
#### Assumptions:
#### 1. baseline CLM: potential NPP
#### 2. Fixed wood NC
#### 3. Even simplier version (assuming 1 soil pool)
####
################################################################################
#### Functions
Perform_Analytical_Run91 <- function(f.flag = 1) {
    #### Function to perform analytical run 91 simulations
    #### eDF: stores equilibrium points
    #### cDF: stores constraint points (curves)
    #### f.flag: = 1 simply plot analytical solution and create individual pdf file
    #### f.flag: = 2 return a list consisting of two dataframes

    ######### Main program
    source("Parameters/Analytical_Run91_Parameters.R")
    
    ### create a range of nc for shoot to initiate
    nfseq <- round(seq(0.001, 0.1, by = 0.0001),5)
    
    ### create nc ratio for wood, root, and allocation coefficients
    a_nf <- as.data.frame(alloc(nfseq))
    
    ### calculate photosynthetic constraint at CO2 = 350 and 700
    P350 <- photo_constraint_full(nf=nfseq, nfdf=a_nf, CO2=CO2_1)

    ### Calculate potential NPP based on VL nutrient recycling constraint
    VL <- VL_constraint_CLM_simplified(a=a_nf)
    
    ### calculate VL equil potential NPP
    VL_act_eq <- solve_VL_CLM_simplified(CO2=CO2_1)

    ### calculate nw and nr for VL equilibrated nf value
    a_act_eq <- alloc(VL_act_eq$nf)
    
    ### calculate soil parameters, e.g. reburial coef.
    s_coef <- soil_coef(df=VL_act_eq$nf, a=a_act_eq)
    
    ### Get omega
    omega_ap <- a_act_eq$af*s_coef$omega_af_pass + a_act_eq$ar*s_coef$omega_ar_pass + a_act_eq$aw*s_coef$omega_aw_pass
    omega_as <- a_act_eq$af*s_coef$omega_af_slow + a_act_eq$ar*s_coef$omega_ar_slow + a_act_eq$aw*s_coef$omega_aw_slow
    
    ### Get C from very-long term nutrient cycling solution
    ### return in g C m-2 
    C_pass_VL_act <- omega_ap*VL_act_eq$NPP/s_coef$decomp_pass/(1-s_coef$qq_pass)*1000.0
    
    ### Calculate long term nutrient constraint
    L <- L_constraint_CLM_simplified(df=nfseq, a=a_nf, 
                                     C_pass=C_pass_VL_act,
                                     Nin_L = Nin)
    
    ### Find long term equilibrium point
    L_act_eq <- solve_L_CLM_simplified(CO2=CO2_1, C_pass=C_pass_VL_act, Nin_L = Nin)
    
    ### allocation stuffs
    a_act_eq <- alloc(L_act_eq$nf)
    
    ### Get Cslow from long nutrient cycling solution
    ### return in g C m-2
    C_slow_L_act <- omega_as*L_act_eq$NPP/s_coef$decomp_slow/(1-s_coef$qq_slow)*1000.0
    
    ### Calculate nutrient release from slow woody pool
    ### return in g N m-2 yr-1
    N_wood_L_act <- a_act_eq$aw*a_act_eq$nw*VL_act_eq$NPP*1000.0
    
    ### Calculate medium term nutrient constraint
    M <- M_constraint_CLM_simplified(df=nfseq,a=a_nf, 
                                     C_pass=C_pass_VL_act, 
                                     C_slow=C_slow_L_act, 
                                     Nin_L = Nin+N_wood_L_act)
    
    ### calculate M equilibrium point
    M_act_eq <- solve_M_CLM_simplified(CO2=CO2_1, 
                                       C_pass=C_pass_VL_act, 
                                       C_slow=C_slow_L_act, 
                                       Nin_L = Nin+N_wood_L_act)
    

    out350DF <- data.frame(CO2_1, nfseq, P350, VL$NPP, 
                           L$NPP, M$NPP)
    colnames(out350DF) <- c("CO2", "nc", "NPP_photo", "NPP_VL",
                            "NPP_L", "NPP_M")
    equil350DF <- data.frame(CO2_1, 
                             VL_act_eq, L_act_eq, M_act_eq)
    colnames(equil350DF) <- c("CO2", 
                              "nc_act_VL", "NPP_act_VL", 
                              "nc_act_L", "NPP_act_L", "nc_act_M", "NPP_act_M")
    
    ##### CO2 = 700
    ### photo constraint
    P700 <- photo_constraint_full(nf=nfseq, nfdf=a_nf, CO2=CO2_2)
    
    ### VL equilibrated point with eCO2
    VL_act_eq <- solve_VL_CLM_simplified(CO2=CO2_2)
    
    ### Find long term equilibrium point
    L_act_eq <- solve_L_CLM_simplified(CO2=CO2_2, C_pass=C_pass_VL_act, Nin_L = Nin)
    
    ### Find medium term equilibrium point
    M_act_eq <- solve_M_CLM_simplified(CO2=CO2_2, 
                                       C_pass=C_pass_VL_act, 
                                       C_slow=C_slow_L_act, 
                                       Nin_L = Nin+N_wood_L_act)
    
    out700DF <- data.frame(CO2_2, nfseq, P700, 
                           VL$NPP, L$NPP, M$NPP)
    colnames(out700DF) <- c("CO2", "nc", "NPP_photo", "NPP_VL",
                            "NPP_L", "NPP_M")
    
    equil700DF <- data.frame(CO2_2, 
                             VL_act_eq, L_act_eq, M_act_eq)
    colnames(equil700DF) <- c("CO2",
                              "nc_act_VL", "NPP_act_VL", 
                              "nc_act_L", "NPP_act_L", "nc_act_M", "NPP_act_M")
 
    ### get the point instantaneous NPP response to doubling of CO2
    df700 <- as.data.frame(cbind(round(nfseq,3), P700))
    inst700 <- inst_NPP(equil350DF$nc_act_VL, df700)
    equil350DF$NPP_I <- inst700$equilNPP
    equil700DF$NPP_I <- inst700$equilNPP
    
    if (f.flag == 1) {
        
        ### plot 2-d plots of nf vs. npp and nf vs. pf
        tiff("Plots/Analytical_Run91_2d.tiff",
             width = 5, height = 5, units = "in", res = 300)
        par(mar=c(5.1,6.1,2.1,2.1))
        
        ### shoot nc vs. NPP
        plot(out350DF$nc, out350DF$NPP_photo, xlim=c(0.001, 0.01),
              ylim=c(0.5, 2.0), 
             type = "l", xlab = "Leaf N:C ratio", 
             ylab = expression(paste("NPP [kg C ", m^-2, " ", yr^-1, "]")),
             col="cyan", lwd = 3, cex.lab=1.5)

        abline(h = seq(0.5, 3.0, 0.5), v = seq(0.01, 0.05, 0.01), col="lightgray", lty = 3)
        
        points(out350DF$nc, out350DF$NPP_VL, type="l", col="tomato", lwd = 3)
        
        points(equil350DF$nc_act_VL, equil350DF$NPP_act_VL, type="p", pch = 19, col = "blue", cex = 2)
        
        points(out350DF$nc, out350DF$NPP_L, type='l',col="violet", lwd = 3)
        
        points(out350DF$nc, out350DF$NPP_M, type="l", col="darkred", lwd = 3)
        
        points(out700DF$nc, out700DF$NPP_photo, col="green", type="l", lwd = 3)
        
        points(equil350DF$nc_act_VL, inst700$equilNPP, type="p", col = "darkgreen", pch=19, cex = 2)
        
        points(equil700DF$nc_act_VL, equil700DF$NPP_act_VL, type="p", col="orange", pch = 19, cex = 2)
        
        points(equil700DF$nc_act_L, equil700DF$NPP_act_L,type="p", col="red", pch = 19, cex = 2)
        
        points(equil700DF$nc_act_M, equil700DF$NPP_act_M, type="p", col="purple", pch = 19, cex = 2)
        
        legend("bottomright", c("P350", "P700", "VL", "L", "M",
                            "A", "B", "C", "D", "E"),
               col=c("cyan","green", "tomato", "violet","darkred","blue", "darkgreen","purple","red", "orange"), 
               lwd=c(2,2,2,2,2,NA,NA,NA,NA,NA),
               pch=c(NA,NA,NA,NA,NA,19,19,19,19,19), cex = 1.0, 
               bg = adjustcolor("grey", 0.8), ncol=2)
        
        dev.off()
        
    } else if (f.flag == 2) {
        
        my.list <- list(cDF = data.frame(rbind(out350DF, out700DF)), 
                        eDF = data.frame(rbind(equil350DF, equil700DF)))
        
        return(my.list)
    } 
}

