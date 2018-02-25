
### Function for nutrient N constraint in medium term 
# not considering exudation
NConsMedium <- function(df, a, Cpass, Cslow, NinL) {
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
    pass <- slow_pool(df, a)
    omegap <- a$af*pass$omegafp + a$ar*pass$omegarp 
    omegas <- a$af*pass$omegafs + a$ar*pass$omegars 
    
    # equation for N constraint with passive, wood, and leaching
    Npass <- (1-pass$qpq) * pass$decomp_p * Cpass * ncp
    Nslow <- pass$decomp_s * Cslow * (1-pass$qsq) * ncs
    
    U0 <- NinL + Npass + Nslow   # will be a constant if decomp rate is constant
    nwood <- a$aw*a$nw
    nburial <- omegap*ncp + omegas*ncs
    nleach <- leachn/(1-leachn) * (a$nfl*a$af + a$nr*(a$ar) + a$nw*a$aw)
    
    NPP_NC <- U0 / (nwood + nburial + nleach)   # will be in g C m-2 yr-1
    NPP <- NPP_NC*10^-3 # returned in kg C m-2 yr-1
    
    df <- data.frame(NPP, nwood,nburial,nleach,a$aw)
    return(df)
}