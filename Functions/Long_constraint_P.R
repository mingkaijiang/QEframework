### Function for nutrient P constraint in longterm ie passive, leaching, wood considered
Long_constraint_P <- function(nf, df, a, Cpass, PinL) {
    # parameters : 
    
    
    # passive pool burial 
    pass <- passive(nf, allocn(nf))
    omegap <- a$af*pass$omegaf + a$ar*pass$omegar 
    
    # equation for P constraint with passive, wood, and leaching
    U0 <- PinL + (1-pass$qq) * pass$decomp * Cpass * pcp   # will be a constant if decomp rate is constant
    pwood <- a$aw*a$pw
    pburial <- omegap*pcp
    pleach <- leachp/(1-leachp-k1) * (a$pfl*a$af + a$pr*a$ar + a$pw*a$aw)
    pocc <- (k3/(k2+k3))*(k1/(1-k1-leachp)) * (a$pfl*a$af + a$pr*a$ar + a$pw*a$aw)
    
    
    NPP_PC <- U0 / (pwood + pburial + pleach + pocc)   # will be in g C m-2 yr-1
    NPP <- NPP_PC*10^-3 # returned in kg C m-2 yr-1
    
    df <- data.frame(NPP, pwood,pburial,pleach, pocc, a$aw)
    return(df)   
}