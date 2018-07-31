
### Function for nutrient N constraint in longterm ie passive, leaching, wood considered
L_constraint_Medlyn_Dewar_no_coupling <- function(a) {
 
    # equation for N constraint with passive, wood, and leaching
    U0 <- Nin
    #nwood <- 0.6 * a$nw
    #nburial <- 0.4 * q * ncp
    nwood <- a$aw * a$nw
    nburial <- (a$af + a$ar) * q * ncp

    # in g C m-2 yr-1
    NPP_NC <- U0 / (nwood + nburial)   
    
    # return in kg C m-2 yr-1
    NPP_N <- NPP_NC*10^-3 
    
    out <- data.frame(NPP_N, nwood, nburial)
    
    colnames(out) <- c("NPP", "nwood", "nburial")
    
    return(out)   
}

L_constraint_Medlyn_Dewar_linear_coupling <- function(a) {
    
    # equation for N constraint with passive, wood, and leaching
    U0 <- Nin
    #nwood <- a$af * scaling_alloc * (a$nw - q * ncp)
    #nburial <- q * ncp
    nwood <- a$aw * a$nw
    nburial <- (a$af + a$ar) * q * ncp
    
    # in g C m-2 yr-1
    NPP_NC <- U0 / (nwood + nburial)   
    
    # return in kg C m-2 yr-1
    NPP_N <- NPP_NC*10^-3 
    
    out <- data.frame(NPP_N, nwood, nburial)
    
    colnames(out) <- c("NPP", "nwood", "nburial")
    
    return(out)   
}