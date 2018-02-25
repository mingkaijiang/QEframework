##### MAIN PROGRAM
mcm_model <- function() {
    nfseq <- seq(0.005,0.05,by=0.001)
    # need allocation fractions here
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
    pass <- passive(nf=VLong$equilnf, a=aequil)
    omegap <- aequil$af*pass$omegaf + aequil$ar*pass$omegar
    CpassVLong <- omegap*VLong$equilNPP/pass$decomp/(1-pass$qq)*1000.0
    NrelwoodVLong <- aequil$aw*aequil$nw*VLong$equilNPP*1000
    
    #now plot long-term constraint with this Cpassive
    NCHUGH <- NConsLong(nf = nfseq,a = a_vec, Cpass=CpassVLong, Nin = 1.0+NrelwoodVLong)
    
    # Solve longterm equilibrium
    equil_long_350 <- solveLong(co2=350, Cpass=CpassVLong, Nin = 1.0+NrelwoodVLong)
    equil_long_700 <- solveLong(co2=700, Cpass=CpassVLong, Nin = 1.0+NrelwoodVLong)
    
    # get the point instantaneous NPP response to doubling of CO2
    df700 <- as.data.frame(cbind(round(nfseq,3), PC700))
    ncref <- round(VLong$equilnf,3)
    
    ## locate the intersect between VL nutrient constraint and CO2 = 700
    VLong700 <- solveVLong(co2=700)
    
    # calculate inst NPP
    inst700 <- inst_NPP(VLong$equilnf, df700)
    
    cDF <- data.frame(nfseq, PC350, PC700, NCHUGH$NPP, NCVLONG$NPP)
    colnames(cDF) <- c("nc", "P350", "P700", "L", "VL")
    eDF <- cbind(VLong, equil_long_700, VLong700, inst700)
    colnames(eDF) <- c("nc_VL", "NPP_VL", "nc_L_700", "NPP_L_700", 
                       "nc_VL_700", "NPP_VL_700",
                       "nc_I", "NPP_I")
    
    return(list(cDF=data.frame(cDF), eDF=data.frame(eDF)))
}
