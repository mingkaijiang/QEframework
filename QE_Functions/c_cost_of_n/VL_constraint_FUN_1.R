VL_constraint_FUN_1 <- function(a) {
    ### Calculating actual NPP based on nutrient recycling
    ### parameters:
    ###            f: fraction of plant enters into soil
    ###            nsoil: soil NC ratio
    
    
    # equation for N constraint with just leaching
    U0 <- Nin
    nleach <- leachn/(1-leachn) * (nfdf$nfl*nfdf$af + nfdf$nr*nfdf$ar + nfdf$nw*nfdf$aw)
    rcn <- 1 / (a$nf*a$af + a$nr*a$ar + a$nw*a$aw)
    
    ### mineral N pool return in kg N m-2 
    Nmin <- Nin *10^-3 / (leachn / (1 - leachn)) 
    
    ### calculate passive N uptake kg N m-2 yr-1
    Npass <- Nmin * (et/sd)
    Npass2 <- pmin(Npass, Nmin)
    
    ###
    cost_active <- c_cost_active(cr=croot)
    
    # Calculate NPP, g C m-2 yr-1
    NPP_NC <- (U0 * (1 + rcn/cost) - rcn * npassive)  / (nleach * (1 + rcn/cost) - 1)
    NPP <- NPP_NC*10^-3     # returned in kg C m-2 yr-1
    
    
    
    
    ### FUN model
    NPP_out <- FUN_model_1(nfdf=a, NPP)
    
    return(NPP_out)   
}