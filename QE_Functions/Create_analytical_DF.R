
#### Create df to store analytical outputs
#### equilDF to store equilibrium points of long and very long term intersection points
#### constraintDF to store G vs. NC and PC over various NC and PC values
####
################################################################################

#### to store all the constraint curves
create_constraint_DF <- function() {
    
    # obtain number of model simulations
    f.list <- list.files(path=paste0(getwd(), "/GDAY/analyses"))
    l <- length(f.list)
    
    # create the matrix
    eDF <- matrix(nrow = l*200, ncol = 13)
    eDF <- as.data.frame(eDF)
    colnames(eDF) <- c("Run", "CO2", "NC", "PC_VL", "PC_350_L",
                       "NPP_350", "NPP_VL", "nleach_VL", "NPP_350_L",
                       "nwood_L", "nburial_L", "nleach_L", "aw")
    
    # add Run information into Run column
    eDF$Run <- rep(c(1:l), each=200)
    
    # add CO2 information into CO2 column
    eDF$CO2 <- rep(c(350, 700), each=100)
    
    return(eDF)
}

#### to store all the constraint curves
create_FACE_constraint_DF <- function() {
    
    # create the matrix
    eDF <- matrix(nrow = 2*182, ncol = 13)
    eDF <- as.data.frame(eDF)
    colnames(eDF) <- c("Run", "CO2", "NC", "PC_VL", "PC_350_L",
                       "NPP_350", "NPP_VL", "nleach_VL", "NPP_350_L",
                       "nwood_L", "nburial_L", "nleach_L", "aw")
    
    # add Run information into Run column
    eDF$Run <- rep(c(1:2), each=182)
    
    # add CO2 information into CO2 column
    eDF$CO2 <- rep(c(350, 700), each=91)
    
    return(eDF)
}

##### to create equilDF to store all the analytical equilibrium points
create_equil_DF <- function() {
    
    # obtain number of model simulations
    f.list <- list.files(path=paste0(getwd(), "/GDAY/analyses"))
    l <- length(f.list)
    
    # create the matrix
    cDF <- matrix(nrow = l*2, ncol = 9)
    cDF <- as.data.frame(cDF)
    colnames(cDF) <- c("Run", "CO2", "NC_VL", "PC_VL","NPP_VL",
                       "NC_L", "PC_L", "NPP_L", "NPP_inst")
    
    # add Run information into Run column
    cDF$Run <- rep(c(1:l), each=2)
    
    # add CO2 information
    cDF$CO2 <- rep(c(350,700), each = 1)
    
    return(cDF)
}

##### to create equilDF to store all the analytical equilibrium points
create_FACE_equil_DF <- function() {
    
    # obtain number of model simulations
    l <- 2
    
    # create the matrix
    cDF <- matrix(nrow = l*2, ncol = 9)
    cDF <- as.data.frame(cDF)
    colnames(cDF) <- c("Run", "CO2", "NC_VL", "PC_VL","NPP_VL",
                       "NC_L", "PC_L", "NPP_L", "NPP_inst")
    
    # add Run information into Run column
    cDF$Run <- rep(c(1:l), each=2)
    
    # add CO2 information
    cDF$CO2 <- rep(c(350,700), each = 1)
    
    return(cDF)
}