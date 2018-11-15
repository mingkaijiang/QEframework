
### Function for nutrient N constraint in longterm ie passive, leaching, wood considered
L_constraint_expl_min <- function(df, a, C_pass, Nin_L) {
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
    pass <- soil_coef(df, a)
    omega_ap <- a$af*pass$omega_af_pass + a$ar*pass$omega_ar_pass + a$aw*pass$omega_aw_pass
    
    # equation for N constraint with passive, wood, and leaching
    U0 <- Nin_L + (1-pass$qq_pass) * pass$decomp_pass * C_pass * ncp   
    nwood <- 0 
    nburial <- omega_ap*ncp
    nleach <- leachn/(1-leachn) * (a$nfl*a$af + a$nr*(a$ar) + a$nw*a$aw)
    
    # in g C m-2 yr-1
    NPP_NC <- (nup * U0) / ((a$nfl*a$af + a$nr*a$ar + a$nw*a$aw) * leachn + nup * nwood + nup * nburial)
    
    # return in kg C m-2 yr-1
    NPP_N <- NPP_NC*10^-3 
    
    # return nleach
    nleach <- leachn * (NPP_NC * (a$nfl*a$af + a$nr*a$ar + a$nw*a$aw)) /nup
    
    out <- data.frame(NPP_N, nwood, nburial, nleach)
    
    colnames(out) <- c("NPP", "nwood", "nburial", "nleach")
    
    return(out)   
}
