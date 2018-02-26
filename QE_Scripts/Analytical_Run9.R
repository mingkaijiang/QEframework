
#### Analytical script Run 9
####
#### Assumptions:
#### 1. baseline CLM: potential NPP
#### 2. Variable wood NC
#### 
####
################################################################################
#### Functions
Perform_Analytical_Run9 <- function(f.flag = 1) {
    #### Function to perform analytical run 9 simulations
    #### eDF: stores equilibrium points
    #### cDF: stores constraint points (curves)
    #### f.flag: = 1 simply plot analytical solution and create individual pdf file
    #### f.flag: = 2 return a list consisting of two dataframes

    ######### Main program
    source("Parameters/Analytical_Run9_Parameters.R")
    
    ### create a range of nc for shoot to initiate
    nfseq <- round(seq(0.001, 0.1, by = 0.001),5)
    
    ### create nc ratio for wood, root, and allocation coefficients
    a_nf <- as.data.frame(alloc(nfseq))
    
    ### calculate photosynthetic constraint at CO2 = 350, which is potential NPP (I think)
    P350 <- photo_constraint_full(nf=nfseq, nfdf=a_nf, CO2=CO2_1)

    ### Calculate potential NPP based on VL nutrient recycling constraint
    VL <- VL_constraint(nf=nfseq, nfdf=a_nf)
    
    ### calculate VL equil potential NPP
    VL_pot_eq <- solve_VL_full(CO2=CO2_1)

    ### calculate nw and nr for VL equilibrated nf value
    a_pot_eq <- alloc(VL_pot_eq$nf)
    
    ### Calculate actual NPP for VL term
    VL_act_eq <- VL_constraint_baseline_CLM_actual(a=a_pot_eq)
    
    ### calculate soil parameters, e.g. reburial coef.
    s_coef <- soil_coef(df=VL_act_eq$nf, a=a_pot_eq)
    
    ### Get omega
    omega_ap <- a_pot_eq$af*s_coef$omega_af_pass + a_pot_eq$ar*s_coef$omega_ar_pass
    omega_as <- a_pot_eq$af*s_coef$omega_af_slow + a_pot_eq$ar*s_coef$omega_ar_slow 
    
    ### Get C from very-long term nutrient cycling solution
    ### return in g C m-2 
    C_pass_VL_pot <- omega_ap*VL_pot_eq$NPP/s_coef$decomp_pass/(1-s_coef$qq_pass)*1000.0
    C_pass_VL_act <- omega_ap*VL_act_eq$NPP/s_coef$decomp_pass/(1-s_coef$qq_pass)*1000.0
    
    ### Calculate long term nutrient constraint
    L <- L_constraint(df=nfseq, a=a_nf, 
                      C_pass=C_pass_VL_pot,
                      Nin_L = Nin)
    
    ### Find long term equilibrium point
    L_pot_eq <- solve_L_full(CO2=CO2_1, C_pass=C_pass_VL_pot, Nin_L = Nin)
    
    ### allocation stuffs
    a_pot_eq <- alloc(L_pot_eq$nf)
    
    ### Find actual L term equilibrium point
    L_act_eq <- L_constraint_baseline_CLM_actual(df=L_pot_eq$nf, a=a_pot_eq, 
                                                 C_pass=C_pass_VL_act,
                                                 Nin_L = Nin)
    
    ### Get Cslow from long nutrient cycling solution
    ### return in g C m-2
    C_slow_L_pot <- omega_as*L_pot_eq$NPP/s_coef$decomp_slow/(1-s_coef$qq_slow)*1000.0
    C_slow_L_act <- omega_as*L_act_eq$NPP/s_coef$decomp_slow/(1-s_coef$qq_slow)*1000.0
    
    ### Calculate nutrient release from slow woody pool
    ### return in g N m-2 yr-1
    N_wood_L_pot <- a_pot_eq$aw*a_pot_eq$nw*VL_pot_eq$NPP*1000.0
    N_wood_L_act <- a_pot_eq$aw*a_pot_eq$nw*VL_act_eq$NPP*1000.0
    
    ### Calculate medium term nutrient constraint
    M <- M_constraint(df=nfseq,a=a_nf, 
                      C_pass=C_pass_VL_pot, 
                      C_slow=C_slow_L_pot, 
                      Nin_L = Nin+N_wood_L_pot)
    
    ### calculate M equilibrium point
    M_pot_eq <- solve_M_full(CO2=CO2_1, 
                             C_pass=C_pass_VL_pot, 
                             C_slow=C_slow_L_pot, 
                             Nin_L = Nin+N_wood_L_pot)
    
    ### allocation stuffs
    a_pot_eq <- alloc(L_pot_eq$nf)
    
    ### actual M NPP equilibrium point
    M_act_eq <- M_constraint_baseline_CLM_actual(df=M_pot_eq$nf,a=a_pot_eq, 
                                                 C_pass=C_pass_VL_act, 
                                                 C_slow=C_slow_L_act, 
                                                 Nin_L = Nin+N_wood_L_act)

    out350DF <- data.frame(CO2_1, nfseq, P350, VL$NPP, 
                           L$NPP, M$NPP)
    colnames(out350DF) <- c("CO2", "nc", "NPP_photo", "NPP_VL",
                            "NPP_L", "NPP_M")
    equil350DF <- data.frame(CO2_1, VL_pot_eq, L_pot_eq, M_pot_eq, 
                             VL_act_eq, L_act_eq, M_act_eq)
    colnames(equil350DF) <- c("CO2", "nc_pot_VL", "NPP_pot_VL", 
                              "nc_pot_L", "NPP_pot_L", "nc_pot_M", "NPP_pot_M",
                              "nc_act_VL", "NPP_act_VL", 
                              "nc_act_L", "NPP_act_L", "nc_act_M", "NPP_act_M")
    
    ##### CO2 = 700
    ### photo constraint
    P700 <- photo_constraint_full(nf=nfseq, nfdf=a_nf, CO2=CO2_2)
    
    ### VL equilibrated point with eCO2
    VL_pot_eq <- solve_VL_full(CO2=CO2_2)
    
    ### calculate nw and nr for VL equilibrated nf value
    a_pot_eq <- alloc(VL_pot_eq$nf)
    
    ### Calculate actual NPP for VL term
    VL_act_eq <- VL_constraint_baseline_CLM_actual(a=a_pot_eq)
    
    ### Find long term equilibrium point
    L_pot_eq <- solve_L_full(CO2=CO2_2, C_pass=C_pass_VL_pot, Nin_L = Nin)
    
    ### Find actual L term equilibrium point
    L_act_eq <- L_constraint_baseline_CLM_actual(df=L_pot_eq$nf, a=a_pot_eq, 
                                                 C_pass=C_pass_VL_act,
                                                 Nin_L = Nin)
    ### Find medium term equilibrium point
    M_pot_eq <- solve_M_full(CO2=CO2_2, 
                             C_pass=C_pass_VL_pot, 
                             C_slow=C_slow_L_pot, 
                             Nin_L = Nin+N_wood_L_pot)
    
    ### Find actual M term point
    M_act_eq <- M_constraint_baseline_CLM_actual(df=M_pot_eq$nf,a=a_pot_eq, 
                                                 C_pass=C_pass_VL_act, 
                                                 C_slow=C_slow_L_act, 
                                                 Nin_L = Nin+N_wood_L_act)
    
    out700DF <- data.frame(CO2_2, nfseq, P700, 
                           VL$NPP, L$NPP, M$NPP)
    colnames(out700DF) <- c("CO2", "nc", "NPP_photo", "NPP_VL",
                            "NPP_L", "NPP_M")
    
    equil700DF <- data.frame(CO2_2, VL_pot_eq, L_pot_eq, M_pot_eq, 
                             VL_act_eq, L_act_eq, M_act_eq)
    colnames(equil700DF) <- c("CO2", "nc_pot_VL", "NPP_pot_VL", 
                              "nc_pot_L", "NPP_pot_L", "nc_pot_M", "NPP_pot_M",
                              "nc_act_VL", "NPP_act_VL", 
                              "nc_act_L", "NPP_act_L", "nc_act_M", "NPP_act_M")
 
    ### get the point instantaneous NPP response to doubling of CO2
    df700 <- as.data.frame(cbind(round(nfseq,3), P700))
    inst700 <- inst_NPP(equil350DF$nc_pot_VL, df700)
    equil350DF$NPP_I <- inst700$equilNPP
    equil700DF$NPP_I <- inst700$equilNPP
    
    if (f.flag == 1) {
        
        ### plot 2-d plots of nf vs. npp and nf vs. pf
        tiff("Plots/Analytical_Run9_2d.tiff",
             width = 5, height = 5, units = "in", res = 300)
        par(mar=c(5.1,6.1,2.1,2.1))
        
        ### shoot nc vs. NPP
        plot(out350DF$nc, out350DF$NPP_photo, xlim=c(0.001, 0.03),
              ylim=c(1.0, 2.0), 
             type = "l", xlab = "Shoot N:C ratio", 
             ylab = expression(paste("Production [kg C ", m^-2, " ", yr^-1, "]")),
             col="cyan", lwd = 3, cex.lab=1.5)

        abline(h = seq(0.5, 3.0, 0.5), v = seq(0.01, 0.05, 0.01), col="lightgray", lty = 3)
        
        points(out350DF$nc, out350DF$NPP_VL, type="l", col="tomato", lwd = 3)
        
        points(equil350DF$nc_pot_VL, equil350DF$NPP_pot_VL, type="p", pch = 21, col = "blue", cex = 2)
        
        points(equil350DF$nc_act_VL, equil350DF$NPP_act_VL, type="p", pch = 19, col = "blue", cex = 2)
        
        points(out350DF$nc, out350DF$NPP_L, type='l',col="violet", lwd = 3)
        
        points(out350DF$nc, out350DF$NPP_M, type="l", col="darkred", lwd = 3)
        
        points(out700DF$nc, out700DF$NPP_photo, col="green", type="l", lwd = 3)
        
        points(equil350DF$nc_pot_VL, inst700$equilNPP, type="p", col = "darkgreen", pch=21, cex = 2)
        
        points(equil700DF$nc_pot_VL, equil700DF$NPP_pot_VL, type="p", col="orange", pch = 21, cex = 2)
        
        points(equil700DF$nc_act_VL, equil700DF$NPP_act_VL, type="p", col="orange", pch = 19, cex = 2)
        
        points(equil700DF$nc_pot_L, equil700DF$NPP_pot_L,type="p", col="red", pch = 21, cex = 2)
        
        points(equil700DF$nc_act_L, equil700DF$NPP_act_L,type="p", col="red", pch = 19, cex = 2)
        
        points(equil700DF$nc_pot_M, equil700DF$NPP_pot_M, type="p", col="purple", pch = 21, cex = 2)
        
        points(equil700DF$nc_act_M, equil700DF$NPP_act_M, type="p", col="purple", pch = 19, cex = 2)
        
        legend("bottomright", c("P350", "P700", "VL", "L", "M",
                            "A1", "B1", "C1", "D1", "E1", "A2", "B2", "C2", "D2", "E2"),
               col=c("cyan","green", "tomato", "violet","darkred","blue", "darkgreen","purple","red", "orange",
                     "blue", "darkgreen","purple","red", "orange"), 
               lwd=c(2,2,2,2,2,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),
               pch=c(NA,NA,NA,NA,NA,21,21,21,21,21,19,19,19,19,19), cex = 1.0, 
               bg = adjustcolor("grey", 0.8), ncol=2)
        
        dev.off()
        
    } else if (f.flag == 2) {
        
        my.list <- list(cDF = data.frame(rbind(out350DF, out700DF)), 
                        eDF = data.frame(rbind(equil350DF, equil700DF)))
        
        return(my.list)
    } 
}

