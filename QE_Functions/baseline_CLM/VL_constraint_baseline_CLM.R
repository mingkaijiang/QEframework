VL_constraint_baseline_CLM_actual <- function(a) {
    ### Calculating actual NPP based on nutrient recycling
    ### parameters:
    ###            f: fraction of plant enters into soil
    ###            nsoil: soil NC ratio

    # equation for N constraint with just leaching
    U0 <- Nin
    
    # Calculate Nmin
    #Nmin <- Nin / (leachn / (1 - leachn))
    Nmin <- Nin / leachn
    
    nplant <- a$nfl*a$af + a$nr*a$ar + a$nw*a$aw
    
    # Calculate NPP
    # NPP_act <- Nmin / (nplant + f * (nsoil - nplant))
    NPP_act <- Nmin / (nplant)
    
    # return in kg C m-2 yr-1
    NPP <- NPP_act *10^-3    
    
    # out df
    out <- data.frame(a$nf, NPP)
    colnames(out) <- c("nf", "NPP")
    
    return(out)   
}