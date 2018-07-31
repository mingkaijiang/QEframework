
#### Analytical script Run 26
####
#### Assumptions:
#### 1. baseline model
#### 2. Fixed
#### 3. baseline N cycle
#### 4. Medlyn and Dewar, 1996, no coupling between allocation leaf and wood
#### 5. simplified soil recycling constraint
####
################################################################################
#### Functions
Perform_Analytical_Run26 <- function(f.flag = 1) {
    #### Function to perform analytical run 26 simulations
    #### eDF: stores equilibrium points
    #### cDF: stores constraint points (curves)
    #### f.flag: = 1 simply plot analytical solution and create individual pdf file
    #### f.flag: = 2 return a list consisting of two dataframes

    ######### Main program
    source("Parameters/Analytical_Run26_Parameters.R")
    
    ### create a range of nc for shoot to initiate
    nfseq <- round(seq(0.001, 0.1, by = 0.001),5)
    
    ### create nc ratio for wood, root, and allocation coefficients
    a_nf <- as.data.frame(alloc_Medlyn_Dewar_no_coupling(nfseq))
    
    ### calculate photosynthetic constraint at CO2 = 350
    P350 <- photo_constraint_full(nf=nfseq, nfdf=a_nf, CO2=CO2_1)

    ### Calculate long term nutrient constraint
    L <- L_constraint_Medlyn_Dewar_no_coupling(a=a_nf)
    
    ### Find long term equilibrium point
    L_eq <- solve_L_Medlyn_Dewar_no_coupling_new(CO2=CO2_1)

    out350DF <- data.frame(CO2_1, nfseq, P350, L$NPP)
    colnames(out350DF) <- c("CO2", "nc", "NPP_photo", "NPP_L")
    equil350DF <- data.frame(CO2_1, L_eq)
    colnames(equil350DF) <- c("CO2", "nc_L", "NPP_L")
    
    ##### CO2 = 700
    ### photo constraint
    P700 <- photo_constraint_full(nf=nfseq, nfdf=a_nf, CO2=CO2_2)
    
    ### Find long term equilibrium point
    L_eq <- solve_L_Medlyn_Dewar_no_coupling_new(CO2=CO2_2)
    
    
    out700DF <- data.frame(CO2_2, nfseq, P700, L$NPP)
    colnames(out700DF) <- c("CO2", "nc", "NPP_photo", "NPP_L")
    
    equil700DF <- data.frame(CO2_2, L_eq)
    colnames(equil700DF) <- c("CO2", "nc_L", "NPP_L")
 
    ### Out put
    if (f.flag == 1) {
  
        
    } else if (f.flag == 2) {
        
        my.list <- list(cDF = data.frame(rbind(out350DF, out700DF)), 
                        eDF = data.frame(rbind(equil350DF, equil700DF)))
        
        return(my.list)
    } 
}

