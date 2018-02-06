
#### Analytical script Run 1
####
#### Assumptions:
#### 1. baseline model
#### 2. Variable wood NC
#### 3. baseline N cycle
####
################################################################################
#### Functions
Perform_Analytical_Run1 <- function(f.flag = 1) {
    #### Function to perform analytical run 1 simulations
    #### eDF: stores equilibrium points
    #### cDF: stores constraint points (curves)
    #### f.flag: = 1 simply plot analytical solution and create individual pdf file
    #### f.flag: = 2 return a list consisting of two dataframes

    ######### Main program
    source("Parameters/Analytical_Run1_Parameters.R")
    
    ### create a range of nc for shoot to initiate
    nfseq <- round(seq(0.001, 0.1, by = 0.001),5)
    
    ### create nc ratio for wood, root, and allocation coefficients
    a_nf <- as.data.frame(alloc(nfseq))
    
    ### calculate photosynthetic constraint at CO2 = 350
    P350 <- photo_constraint_full(nf=nfseq, nfdf=a_nf, CO2=CO2_1)
    P700 <- photo_constraint_full(nf=nfseq, nfdf=a_nf, CO2=CO2_2)
    
    ### calculate very long term NC constraint on NPP, respectively
    VL <- VL_constraint(nf=nfseq, nfdf=a_nf)
    
    ### finding the equilibrium point between photosynthesis and very long term nutrient constraints
    VL_eq <- solve_VL_full(CO2=CO2_1)
    
    ### calculate nw and nr for VL equilibrated nf value
    a_eq <- alloc(VL_eq$nf)
    
    ### calculate soil parameters, e.g. reburial coef.
    s_coef <- soil_coef(df=VL_eq$nf, a=a_eq)
    omega_ap <- a_eq$af*s_coef$omega_af_pass + a_eq$ar*s_coef$omega_ar_pass
    omega_as <- a_eq$af*s_coef$omega_af_slow + a_eq$ar*s_coef$omega_ar_slow
    
    ### Get C from very-long term nutrient cycling solution
    ### return in g C m-2 
    C_pass_VL <- omega_ap*VL_eq$NPP/s_coef$decomp_pass/(1-s_coef$qq_pass)*1000.0

    ### Calculate long term nutrient constraint
    L <- L_constraint(df=nfseq, a=a_nf, 
                      C_pass=C_pass_VL,
                      Nin_L = Nin)
    
    ### Find long term equilibrium point
    L_eq <- solve_L_full(CO2=CO2_1, C_pass=C_pass_VL, Nin_L = Nin)
    
    ### Get Cslow from long nutrient cycling solution
    ### return in g C m-2
    C_slow_L <- omega_as*L_eq$NPP/s_coef$decomp_slow/(1-s_coef$qq_slow)*1000.0
    
    ### Calculate nutrient release from slow woody pool
    NrelwoodVLong <- aequiln$aw*aequiln$nw*VLong_equil$equilNPP*1000.0
    
    ### Calculate medium term nutrient constraint
    NCMEDIUM <- NConsMedium(df=nfseq, 
                            a=a_nf, 
                            Cpass=CpassVLong, 
                            Cslow=CslowLong, 
                            NinL = Nin+NrelwoodVLong)
    
    Medium_equil_350 <- solveMedium(CO2=CO2_1, 
                                    Cpass=CpassVLong, 
                                    Cslow=CslowLong, 
                                    NinL = Nin+NrelwoodVLong)
    

    out350DF <- data.frame(CO2_1, nfseq, Photo350, NCVLONG$NPP_N, 
                           NCLONG$NPP, NCMEDIUM$NPP)
    colnames(out350DF) <- c("CO2", "nc", "NPP_photo", "NPP_VL",
                            "NPP_L", "NPP_M")
    equil350DF <- data.frame(CO2_1, VLong_equil, Long_equil, Medium_equil_350)
    colnames(equil350DF) <- c("CO2", "nc_VL", "NPP_VL", 
                              "nc_L", "NPP_L", "nc_M", "NPP_M")
    
    ##### CO2 = 700
    ### N:C ratio
    nfseq <- round(seq(0.001, 0.1, by = 0.001),5)
    a_nf <- as.data.frame(allocn(nfseq))
    
    ### calculate NC vs. NPP at CO2 = 350 respectively
    Photo700 <- photo_constraint_full_cn(nfseq, a_nf, CO2_2)
    
    ### calculate very long term NC and PC constraint on NPP, respectively
    NCVLONG <- VLong_constraint_N(nf=nfseq, nfdf=a_nf)
    
    ### finding the equilibrium point between photosynthesis and very long term nutrient constraints
    VLong_equil <- solveVLong_full_cn(CO2=CO2_2)
    
    ### Find long term equilibrium point
    Long_equil <- solveLong_full_cn(CO2=CO2_2, Cpass=CpassVLong, NinL = Nin)
    
    ### Find medium term equilibrium point
    Medium_equil_700 <- solveMedium(CO2_2, Cpass = CpassVLong, Cslow = CslowLong, 
                                             NinL=Nin+NrelwoodVLong)
    
    out700DF <- data.frame(CO2_2, nfseq, Photo700, 
                           NCVLONG$NPP_N, NCLONG$NPP, NCMEDIUM$NPP)
    colnames(out700DF) <- c("CO2", "nc", "NPP_photo", "NPP_VL",
                            "NPP_L", "NPP_M")
    
    equil700DF <- data.frame(CO2_2, VLong_equil, Long_equil, Medium_equil_700)
    colnames(equil700DF) <- c("CO2", "nc_VL", "NPP_VL", 
                              "nc_L", "NPP_L", "nc_M", "NPP_M")
 
    ### get the point instantaneous NPP response to doubling of CO2
    df700 <- as.data.frame(cbind(round(nfseq,3), Photo700))
    inst700 <- inst_NPP(equil350DF$nc_VL, df700)
    
    if (f.flag == 1) {
        
        ### plot 2-d plots of nf vs. npp and nf vs. pf
        tiff("Plots/Analytical_Run1_2d.tiff",
             width = 5, height = 5, units = "in", res = 300)
        par(mar=c(5.1,6.1,2.1,2.1))
        
        ### nice plot
        # require(plotly)
        # plot_ly(data = out350DF, x = ~nc, y = ~NPP_photo,
        #              name = "P350", type = "scatter", mode = "lines") %>%
        #     add_trace(y = ~NPP_VL, name = "VL", mode = "lines") %>%
        #     add_trace(y = ~NPP_L, name = "L", mode = "lines") %>%
        #     add_trace(y = ~NPP_M, name = "M", mode = "lines") %>%
        #     layout(xaxis = list(range = c(0.01, 0.05)),
        #            yaxis = list(range = c(0, 3)))
        
        ### shoot nc vs. NPP
        plot(out350DF$nc, out350DF$NPP_photo, xlim=c(0.018, 0.030),
              ylim=c(1.605, 1.7), 
             type = "l", xlab = "Shoot N:C ratio", 
             ylab = expression(paste("Production [kg C ", m^-2, " ", yr^-1, "]")),
             col="cyan", lwd = 3, cex.lab=1.5)
        
        points(equil350DF$nc_VL, equil350DF$NPP_L, col="black", lty = 3)
        points(equil350DF$nc_VL, equil350DF$NPP_M, col="red", lty = 3)
        points(equil700DF$nc_VL, equil700DF$NPP_M, col="red", lty = 3)
        
        abline(h = seq(0.5, 3.0, 0.5), v = seq(0.01, 0.05, 0.01), col="lightgray", lty = 3)
        
        points(out350DF$nc, out350DF$NPP_VL, type="l", col="tomato", lwd = 3)
        
        points(equil350DF$nc_VL, equil350DF$NPP_VL, type="p", pch = 19, col = "blue", cex = 2)
        
        points(out350DF$nc, out350DF$NPP_L, type='l',col="violet", lwd = 3)
        
        points(nfseq, NCMEDIUM$NPP, type="l", col="darkred", lwd = 3)
        
        points(out700DF$nc, out700DF$NPP_photo, col="green", type="l", lwd = 3)
        
        points(equil350DF$nc_VL, inst700$equilNPP, type="p", col = "darkgreen", pch=19, cex = 2)
        
        points(equil700DF$nc_VL, equil700DF$NPP_VL, type="p", col="orange", pch = 19, cex = 2)
        
        points(equil700DF$nc_L, equil700DF$NPP_L,type="p", col="red", pch = 19, cex = 2)
        
        points(Medium_equil_700$equilnf, Medium_equil_700$equilNPP, type="p", col="purple", pch = 19, cex = 2)
        
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

