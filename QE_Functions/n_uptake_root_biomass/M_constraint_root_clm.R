### Function for nutrient N constraint in medium term ie passive, leaching, wood considered
# specifically for uptake as a function of root  - O-CN approach
# i.e. N uptake as a saturating function of mineral N
M_constraint_root_clm <- function(df, a, C_pass, C_slow, Nin_L) {
    
    # nc scalar
    scalar_n <- (1/df - cn.leaf.min) / (cn.leaf.max - cn.leaf.min)
    
    # passive pool burial 
    s_coef <- soil_coef(df, a)
    omega_ap <- a$af*s_coef$omega_af_pass + a$ar*s_coef$omega_ar_pass + a$aw*s_coef$omega_aw_pass
    omega_as <- a$af*s_coef$omega_af_slow + a$ar*s_coef$omega_ar_slow + a$aw*s_coef$omega_aw_slow
    
    # N comes from decomposition of pass and slow pools
    Npass <- (1-s_coef$qq_pass) * s_coef$decomp_pass * C_pass * ncp
    Nslow <- (1-s_coef$qq_slow) * s_coef$decomp_slow * C_slow * ncs
    
    # equation for N constraint with passive, wood, and leaching
    U0 <- Nin_L + Npass + Nslow   
    nwood <- a$aw*a$nw
    nburial <- omega_ap*ncp + omega_as*ncs
 
    # n leach
    nleach <- leachn * ksmin * (a$af * a$nf + a$aw * a$nw + a$ar * 
                                a$nr) / (umax * (a$ar/sr) * scalar_temp * scalar_n - (a$af * a$nf + a$aw * a$nw + a$ar * a$nr))
    
    #browser()
    
    # in g C m-2 yr-1
    NPP_NC <- (U0 - nleach) / (nwood + nburial)
    
    # returned in kg C m-2 yr-1
    NPP <- NPP_NC*10^-3 
    
    out <- data.frame(NPP, nwood, nburial, nleach)
    
    return(out)   
}