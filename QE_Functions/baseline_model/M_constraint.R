
### Function for nutrient N constraint in medium term 
# not considering exudation
M_constraint <- function(df, a, C_pass, C_slow, Nin_L) {
    # passed are df and a, the allocation and plant N:C ratios
    # parameters : 
    # Nin is fixed N inputs (N deposition annd fixation) in g m-2 yr-1 (could vary fixation)
    # nleach is the rate of n leaching of the mineral pool (per year)
    # Tsoil is effective soil temperature for decomposition
    # Texture is the fine soil fraction
    # ligfl and ligrl are the lignin:C fractions in the foliage and root litter
    # Cpass is the passive pool size in g C m-2
    # ncp is the NC ratio of the passive pool in g N g-1 C
    
    # passive pool burial 
    s_coef <- soil_coef(df, a)
    omega_ap <- a$af*s_coef$omega_af_pass + a$ar*s_coef$omega_ar_pass 
    omega_as <- a$af*s_coef$omega_af_slow + a$ar*s_coef$omega_ar_slow 
    
    # equation for N constraint with passive, wood, and leaching
    Npass <- (1-s_coef$qq_pass) * s_coef$decomp_pass * C_pass * ncp
    Nslow <- (1-s_coef$qq_slow) * s_coef$decomp_slow * C_slow * ncs
    
    U0 <- Nin_L + Npass + Nslow   
    nwood <- a$aw*a$nw
    nburial <- omega_ap*ncp + omega_as*ncs
    nleach <- leachn/(1-leachn) * (a$nfl*a$af + a$nr*(a$ar) + a$nw*a$aw)
    
    ### return g C m-2 yr-1
    NPP_NC <- U0 / (nwood + nburial + nleach)   
    
    ### return kg C m-2 yr-1
    NPP <- NPP_NC*10^-3 
    
    df <- data.frame(NPP, nwood, nburial, nleach)
    return(df)
}
