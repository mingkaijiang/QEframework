
L_constraint_root_clm <- function(df, a, C_pass, Nin_L) {
    ### Function for nutrient N constraint in longterm ie passive, leaching, wood considered
    # specifically for uptake as a function of root  - O-CN approach
    # i.e. N uptake as a saturating function of mineral N
    # passed are df and a, the allocation and plant N:C ratios
    # parameters : 
    # Nin is fixed N inputs (N deposition annd fixation) in g m-2 yr-1 (could vary fixation)
    # nleach is the rate of n leaching of the mineral pool (per year)
    # Tsoil is effective soil temperature for decomposition
    # Texture is the fine soil fraction
    # ligfl and ligrl are the lignin:C fractions in the foliage and root litter
    # Cpass is the passive pool size in g C m-2
    # ncp is the NC ratio of the passive pool in g N g-1 C
    
    # nc scalar
    scalar_n <- (1/df - cn.leaf.min) / (cn.leaf.max - cn.leaf.min)
    
    # passive pool burial 
    pass <- soil_coef(df, a)
    omega_ap <- a$af*pass$omega_af_pass + a$ar*pass$omega_ar_pass 
    
    U0 <- Nin_L + (1-pass$qq_pass) * pass$decomp_pass * C_pass * ncp   
    nburial <- omega_ap*ncp
    nleach <- leachn * ksmin * (a$af * a$nf + a$aw * a$nw + a$ar * 
                                a$nr) / (umax * (a$ar/sr) * scalar_temp * scalar_n
                                         - (a$af * a$nf + a$aw * a$nw + a$ar * a$nr))
    
    #browser()
    
    # will be in g C m-2 yr-1
    NPP_NC <- (U0 - nleach) / (nburial) 
    
    # returned in kg C m-2 yr-1
    NPP <- NPP_NC*10^-3 
    
    out <- data.frame(NPP, nburial, nleach)
    
    return(out)   
}


