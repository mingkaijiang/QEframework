#### Functions
P_limitation_effect <- function() {
    #### Perform CNP only analysis of VL pools
    source("Parameters/Analytical_Run1_Parameters.R")
    
    # create a range of nc for shoot to initiate
    nfseq <- round(seq(0.01, 0.1, by = 0.001),5)
    a_nf <- as.data.frame(allocn(nfseq))
    
    # using very long term relationship to calculate pf from nf
    pfseq <- inferpfVL(nfseq, a_nf)
    a_pf <- as.data.frame(allocp(pfseq))
    
    # calculate photosynthetic constraint at CO2 = 350
    photo_350_cnp <- photo_constraint_full_cnp(nfseq, pfseq, a_nf, a_pf, CO2_1)
    photo_700_cnp <- photo_constraint_full_cnp(nfseq, pfseq, a_nf, a_pf, CO2_2)
    
    ### calculate very long term NC and PC constraint on NPP, respectively
    vlong_cnp <- VLong_constraint_N(nf=nfseq, nfdf=a_nf)
    
    ### finding the equilibrium point between photosynthesis and very long term nutrient constraints
    VLong_equil_cnp <- solveVLong_full_cnp(CO2=CO2_1)
    VLong_equil_cnp_new <- solveVLong_full_cnp(CO2=CO2_2)
    
    #### Perform CN only analysis of VL pools
    source("Parameters/Analytical_Run2_Parameters.R")
    
    # N:C ratios for x-axis
    nfseq <- seq(0.01,0.1,by=0.001)
    # need allocation fractions here
    a_vec <- allocn(nfseq)
    
    # plot photosynthetic constraints
    photo_350_cn <- photo_constraint_full_cn(nfseq,a_vec,CO2=CO2_1)
    photo_700_cn <- photo_constraint_full_cn(nfseq,a_vec,CO2=CO2_2)
    
    #plot very-long nutrient cycling constraint
    vlong_cn <- VLong_constraint_N(nfseq,a_vec)
    
    #solve very-long nutrient cycling constraint
    VLong_equil_cn <- solveVLong_full_cn(CO2=CO2_1)
    VLong_equil_cn_new <- solveVLong_full_cn(CO2=CO2_2)
    
    co2_effect_cnp <- (VLong_equil_cnp_new$equilNPP - VLong_equil_cnp$equilNPP) / VLong_equil_cnp$equilNPP * 100
    co2_effect_cn <- (VLong_equil_cn_new$equilNPP - VLong_equil_cn$equilNPP) / VLong_equil_cn$equilNPP * 100
    
    #### Plotting
    tiff("Plots/Effect_of_P_limitation_on_CO2_fertilization.tiff",
         width = 8, height = 7, units = "in", res = 300)
    par(mar=c(5.1,6.1,2.1,2.1))
    
    # shoot nc vs. NPP
    plot(nfseq, photo_350_cnp, xlim=c(0.01, 0.05),
         ylim=c(0.0, 3.0), 
         type = "l", xlab = "Shoot N:C ratio", 
         ylab = expression(paste("Production [kg C ", m^-2, " ", yr^-1, "]")),
         col="cyan", lwd = 3, cex.lab = 1.5)
    points(nfseq, vlong_cnp$NPP_N, type="l", col="tomato", lwd = 3)
    points(VLong_equil_cnp$equilnf, VLong_equil_cnp$equilNPP, type="p", pch = 19, col = "blue", cex = 2.0)
    points(nfseq, photo_350_cn, type="l", col = "cyan", lty = 3, lwd = 3)
    points(VLong_equil_cn$equilnf, VLong_equil_cn$equilNPP, type="p", pch = 1, col = "blue", cex = 2.0)
    
    points(nfseq, photo_700_cnp, type="l", col="green", lwd=3)
    points(nfseq, photo_700_cn, type="l", col="green", lwd=3, lty=3)
    points(VLong_equil_cnp_new$equilnf, VLong_equil_cnp_new$equilNPP, type="p", pch = 19, col = "orange", cex=2.0)
    points(VLong_equil_cn_new$equilnf, VLong_equil_cn_new$equilNPP, type="p", pch = 1, col = "orange", cex=2.0)
    
    legend("bottomleft", c(expression("NP"[a]), 
                         expression("N"[a]), 
                         expression("NP"[e]), 
                         expression("N"[e]), 
                         "VL", 
                         "A", "B", "C", "D"),
           col=c("cyan","cyan", "green", "green", "tomato",
                 "blue","orange", "blue", "orange"), 
           lwd=c(2,2,2,2,2, NA, NA, NA, NA),  lty=c(1,3,1,3,1, NA, NA, NA, NA),
           pch = c(NA, NA, NA, NA, NA, 19,19,1,1), cex = 1.2, ncol=2)
    
    dev.off()
    
}


#### Program
P_limitation_effect()
