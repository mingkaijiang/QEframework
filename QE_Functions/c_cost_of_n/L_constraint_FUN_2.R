
### Function for nutrient N constraint in longterm ie passive, leaching considered
L_constraint_FUN_2 <- function(df, a, C_pass, Nin_L) {
    
    # passive pool burial 
    pass <- soil_coef(df, a)
    omega_ap <- a$af*pass$omega_af_pass + a$ar*pass$omega_ar_pass 
    
    # equation for N constraint with passive pool
    U0 <- Nin_L + (1-pass$qq_pass) * pass$decomp_pass * C_pass * ncp  
    nburial <- omega_ap*ncp

    nplant <- a$nf*a$af + a$nr*a$ar + a$nw*a$aw
    nleach <- (leachn/(1-leachn)) * (nplant + f * (nsoil - nplant))
        
    # in kg C m-2 yr-1 
    NPP <- U0 * 10^-3 / (nburial + nleach)
    
    # add FUN modifier
    NPP_out <- FUN_model_2(nfdf=a, NPP)
    
    return(NPP_out)   
}
