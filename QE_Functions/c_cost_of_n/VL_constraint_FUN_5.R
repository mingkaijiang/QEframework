VL_constraint_FUN_5 <- function(a) {
    ### Calculating actual NPP based on nutrient recycling
    ### parameters:
    ###            f: fraction of plant enters into soil
    ###            nsoil: soil NC ratio
    
    # Calculate Nmin, g m-2 yr-1
    Nmin <- Nin * 10^-3 / (leachn / (1 - leachn)) 
    
    nplant <- a$nf*a$af + a$nr*a$ar + a$nw*a$aw
    
    # Calculate NPP
    NPP <- Nmin / (nplant + f * (nsoil - nplant))
    
    ### FUN model
    NPP_out <- FUN_model_5(nfdf=a, NPP)
    
    return(NPP_out)   
}