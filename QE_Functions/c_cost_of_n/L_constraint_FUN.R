
### Function for nutrient N constraint in longterm ie passive, leaching considered
L_constraint_FUN <- function(df, a, C_pass, Nin_L) {
    
    # passive pool burial 
    pass <- soil_coef(df, a)
    omega_ap <- a$af*pass$omega_af_pass + a$ar*pass$omega_ar_pass 
    
    # equation for N constraint with passive pool
    U0 <- Nin_L + (1-pass$qq_pass) * pass$decomp_pass * C_pass * ncp   
    nburial <- omega_ap*ncp

    nplant <- a$nfl*a$af + a$nr*a$ar + a$nw*a$aw
    nleach <- (leachn/(1-leachn)) * (nplant + f * (nsoil - nplant))
        
    # in g C m-2 yr-1 
    # NPP_NC <- U0 / ((leachn + 1) * (nplant + f * (nsoil - nplant)))
    NPP_NC <- U0 / (nburial + nleach)

    # return in kg C m-2 yr-1
    NPP <- NPP_NC*10^-3 
    
    # add FUN modifier
    NPP_out <- FUN_model(nfdf=a, NPP)
    
    return(NPP_out)   
}
