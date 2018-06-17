
M_constraint_FUN_6 <- function(df, a, C_pass, C_slow, Nin_L) {
 
    # passive pool burial 
    s_coef <- soil_coef(df, a)
    omega_ap <- a$af*s_coef$omega_af_pass + a$ar*s_coef$omega_ar_pass 
    omega_as <- a$af*s_coef$omega_af_slow + a$ar*s_coef$omega_ar_slow 
    
    # N comes from decomposition of pass and slow pools
    Npass <- (1-s_coef$qq_pass) * s_coef$decomp_pass * C_pass * ncp
    Nslow <- (1-s_coef$qq_slow) * s_coef$decomp_slow * C_slow * ncs
    
    # equation for N constraint with passive, wood, and leaching
    U0 <- Nin_L + Npass + Nslow   
    nwood <- a$aw*a$nw
    nburial <- omega_ap*ncp + omega_as*ncs
    nplant <- a$nf*a$af + a$nr*a$ar + a$nw*a$aw
    nleach <- (leachn/(1-leachn)) * (nplant + f * (nsoil - nplant))
    
    ### return kg C m-2 yr-1
    NPP <- U0 * 10^-3 / (nburial + nleach + nwood)
  
    NPP_out <- data.frame(a$nf, NPP)
    colnames(NPP_out) <- c("nf", "NPP")
    
    return(NPP_out)   
}
