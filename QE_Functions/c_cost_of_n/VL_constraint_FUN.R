VL_constraint_FUN <- function(a) {
    ### Calculating actual NPP based on nutrient recycling
    ### parameters:
    ###            f: fraction of plant enters into soil
    ###            nsoil: soil NC ratio

    # equation for N constraint with just leaching
    U0 <- Nin
    
    # Calculate Nmin
    Nmin <- Nin / (leachn / (1 - leachn))
    # Nmin <- Nin / leachn
    
    nplant <- a$nf*a$af + a$nr*a$ar + a$nw*a$aw
    
    # Calculate NPP
    NPP_act <- Nmin / (nplant + f * (nsoil - nplant))
    
    # return in kg C m-2 yr-1
    NPP <- NPP_act *10^-3   
    
    ### FUN model
    NPP_out <- FUN_model(nfdf=a, NPP)
    
    return(NPP_out)   
}