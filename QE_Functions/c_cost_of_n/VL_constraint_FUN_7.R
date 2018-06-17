VL_constraint_FUN_7 <- function(a) {
    ### Calculating actual NPP based on nutrient recycling
    ### parameters:
    ###            f: fraction of plant enters into soil
    ###            nsoil: soil NC ratio
    
    # Calculate Nmin, g m-2 yr-1
    U0 <- Nin * 10^-3
    Nmin <- U0 / (leachn / (1 - leachn)) 
    
    nplant <- a$nf*a$af + a$nr*a$ar + a$nw*a$aw
    
    # Calculate NPP
    NPP <- Nmin / (nplant + f * (nsoil - nplant))
    
    ### FUN model
    NPP_out <- FUN_model_7(u0=U0, nfdf=a, NPP)
    
    return(NPP_out)   
}