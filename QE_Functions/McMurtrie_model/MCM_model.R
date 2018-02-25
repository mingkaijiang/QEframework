##### MAIN PROGRAM
mcm_model <- function() {
    
    source("Parameters/Analytical_MCM_Parameters.R")
    
    # nc array
    nfseq <- seq(0.01,0.05,by=0.001)
    
    # allocation
    a_vec <- alloc_mcm(nfseq)
    
    # plot photosynthetic constraints
    PC350 <- solvePC(nfseq,a_vec$af,co2=350)
    PC700 <- solvePC(nfseq,a_vec$af,co2=700)

    #plot very-long nutrient cycling constraint
    NCVLONG <- NConsVLong(nf=nfseq,a=a_vec,Nin=1.0)
    
    #solve very-long nutrient cycling constraint
    VLong <- solveVLong(co2=350)
    
    #get Cpassive from very-long nutrient cycling solution
    aequil <- alloc_mcm(VLong$equilnf)
    
    pass <- slow_pool(df=VLong$equilnf, a=aequil)
    omegap <- aequil$af*pass$omegafp + aequil$ar*pass$omegarp
    omegas <- aequil$af*pass$omegafs + aequil$ar*pass$omegars
    
    CpassVLong <- omegap*VLong$equilNPP/pass$decomp_p/(1-pass$qpq)*1000.0
    
    NrelwoodVLong <- aequil$aw*aequil$nw*VLong$equilNPP*1000
    
    #now plot long-term constraint with this Cpassive
    NCHUGH <- NConsLong(df = nfseq,a = a_vec, Cpass=CpassVLong, NinL = 1.0)
    
    # Solve longterm equilibrium
    equil_long_350 <- solveLong(co2=350, Cpass=CpassVLong, NinL = 1.0)
    equil_long_700 <- solveLong(co2=700, Cpass=CpassVLong, NinL = 1.0)
    
    # get the point instantaneous NPP response to doubling of CO2
    df700 <- as.data.frame(cbind(round(nfseq,3), PC700))

    # slow soil pool
    CslowLong <- omegas*equil_long_350$equilNPP/pass$decomp_s/(1-pass$qsq)*1000.0
    
    # plot medium nutrient cycling constraint
    NCMEDIUM <- NConsMedium(nfseq, a_vec, CpassVLong, CslowLong, NinL=1.0+NrelwoodVLong)
    
    # solve medium term equilibrium at CO2 = 700 ppm
    equil_medium_700 <- solveMedium(CO2=700,Cpass=CpassVLong,Cslow=CslowLong,NinL=1.0+NrelwoodVLong)
    
    ## locate the intersect between VL nutrient constraint and CO2 = 700
    VLong700 <- solveVLong(co2=700)
    
    # calculate inst NPP
    inst700 <- inst_NPP(VLong$equilnf, df700)
    
    cDF <- data.frame(nfseq, PC350, PC700, NCHUGH$NPP, NCVLONG$NPP, NCMEDIUM$NPP)
    colnames(cDF) <- c("nc", "P350", "P700", "L", "VL", "M")
    eDF <- cbind(VLong, equil_medium_700, equil_long_700, VLong700, inst700)
    colnames(eDF) <- c("nc_VL", "NPP_VL", "nc_M_700", "NPP_M_700",
                       "nc_L_700", "NPP_L_700", 
                       "nc_VL_700", "NPP_VL_700",
                       "nc_I", "NPP_I")
    
    return(list(cDF=data.frame(cDF), eDF=data.frame(eDF)))
}
