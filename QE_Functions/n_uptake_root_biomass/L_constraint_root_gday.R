
L_constraint_root_gday <- function(df, a, Nin_L, 
                                   C_pass) {
    ### Function for nutrient N constraint in longterm ie passive, leaching, wood considered
    # specifically for uptake as a function of root - GDAY approach
    # i.e. N uptake as a saturating function of root biomass
    # passed are df and a, the allocation and plant N:C ratios
    # parameters : 
    # Nin is fixed N inputs (N deposition annd fixation) in g m-2 yr-1 (could vary fixation)
    # nleach is the rate of n leaching of the mineral pool (per year)
    # Tsoil is effective soil temperature for decomposition
    # Texture is the fine soil fraction
    # ligfl and ligrl are the lignin:C fractions in the foliage and root litter
    # Cpass is the passive pool size in g C m-2
    # ncp is the NC ratio of the passive pool in g N g-1 C
    # sr is the decay rate of root in yr-1
    # kr is the value of root carbon at which 50% of the available N is taken up
    
    # passive pool burial 
    pass <- soil_coef(df, a)
    omega_ap <- a$af*pass$omega_af_pass + a$ar*pass$omega_ar_pass + a$aw*pass$omega_aw_pass
    
    # equation for N constraint with passive, wood, and leaching
    U0 <- Nin_L + (1-pass$qq_pass) * pass$decomp_pass * C_pass * ncp   
    nwood <- 0 
    nburial <- omega_ap*ncp
 
    # allocated N to plant   
    a_plant <- a$nfl*a$af + a$nr*a$ar + a$nw*a$aw
    
    # equation for NPP
    NPP_NC <- (U0 - a_plant * leachn * kr * (sr / a$ar)) / (nwood + nburial + a_plant*leachn)
    
    # returned in kg C m-2 yr-1
    NPP <- NPP_NC*10^-3 
    
    # leaching
    nleach <- leachn * (U0 - nwood * NPP - nburial * NPP)
    
    # out df
    out <- data.frame(NPP, nwood, nburial, nleach)
    
    return(out)   
}

