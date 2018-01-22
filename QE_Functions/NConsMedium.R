
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

### Function for nutrient N constraint in medium term 
# considering priming
NConsMedium_priming <- function(df, a, Cpass, Cslow, NinL) {
    # passed are df and a, the allocation and plant N:C ratios
    # parameters : 
    # Nin is fixed N inputs (N deposition annd fixation) in g m-2 yr-1 (could vary fixation)
    # nleach is the rate of n leaching of the mineral pool (per year)
    # Tsoil is effective soil temperature for decomposition
    # Texture is the fine soil fraction
    # ligfl and ligrl are the lignin:C fractions in the foliage and root litter
    # Cpass is the passive pool size in g C m-2
    # ncp is the NC ratio of the passive pool in g N g-1 C
    
    len <- length(df)
    
    ans <- c()
    
    # burial fractions
    pas <- 0.996 - (0.85-0.68*Texture)
    psa <- 0.42
    ppa <- 0.45
    pap <- 0.004
    psp <- 0.03
    
    
    for (i in 1:len) {
        fPC <- function(NPP) {
            # passive and slow pool burial 
            pass <- slow_priming(df[i], a[i,], NPP)
            omegap <- a[i,]$af*pass$omegafp + a[i,]$ar*pass$omegarp 
            omegas <- a[i,]$af*pass$omegafs + a[i,]$ar*pass$omegars 
            
            # equation for N constraint with passive, slow, wood, and leaching
            Npass <- (1-pass$qpq) * pass$decomp_p * Cpass * ncp
            Nslow <- pass$decomp_s * Cslow * (1-pass$qsq) * ncs

            U0 <- NinL + Npass + Nslow
            nwood <- a[i,]$aw*a[i,]$nw
            nburial <- omegap*ncp + omegas*ncs
            nleach <- leachn / (1-leachn) * (a[i,]$nfl*a[i,]$af + a[i,]$nr*a[i,]$ar 
                                             + a[i,]$nw*a[i,]$aw)
            
            NPP_NC <- U0 / (nleach + nwood + nburial)
            
            out <- NPP_NC*10^-3 - NPP # returned in kg C m-2 yr-1
            
        }
        
        ans[i] <- uniroot(fPC,interval=c(0.1,20), trace=T)$root
        
    }
    
    out <- data.frame(ans, a$ariz)
    colnames(out) <- c("NPP", "ariz")
    return(out)   
}



### Function for nutrient N constraint in medium term ie passive, leaching, wood considered
# specifically for explicit mineral N
NConsMedium_expl_min <- function(df, a, Cpass, Cslow, NinL) {
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
    
    #NPP_NC <- ((U0 - nwood - nburial) / leachn) * nuptakerate / (a$nfl*a$af + a$nr*a$ar + a$nw*a$aw)   # will be in g C m-2 yr-1
    NPP_NC <- (nuptakerate * U0) / ((a$nfl*a$af + a$nr*a$ar + a$nw*a$aw) * leachn + nuptakerate * nwood + nuptakerate * nburial)
    
    NPP <- NPP_NC*10^-3 # returned in kg C m-2 yr-1
    
    nleach <- leachn * (NPP_NC * (a$nfl*a$af + a$nr*a$ar + a$nw*a$aw)) /nuptakerate
    
    df <- data.frame(NPP, nwood,nburial,nleach,a$aw)
    return(df)   
}

### Function for nutrient N constraint in longterm ie passive, leaching, wood considered
# specifically for uptake as a function of root - GDAY approach
# i.e. N uptake as a saturating function of root biomass
NConsMedium_root_gday <- function(df, a, NinL, 
                                Cpass, Cslow) {
    
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
    
    A_NF <- a$nfl*a$af + a$nr*a$ar + a$nw*a$aw
    
    # equation for NPP
    NPP_NC <- (U0 - A_NF * leachn * kr * (sr / a$ar)) / (nwood + nburial + A_NF*leachn)
    
    NPP_N <- NPP_NC*10^-3 # returned in kg C m-2 yr-1
    
    nleach <- leachn * (U0 - nwood * NPP_N - nburial * NPP_N)
    
    df <- data.frame(NPP_N, nwood,nburial,nleach,a$aw)
    return(df)   
}

### Function for nutrient N constraint in medium term ie passive, leaching, wood considered
# specifically for uptake as a function of root  - O-CN approach
# i.e. N uptake as a saturating function of mineral N
NConsMedium_root_ocn <- function(df, a, Cpass, Cslow, NinL) {
    
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
    nleach <- leachn * k * (a$af * a$nf + a$aw * a$nw + a$ar * a$nr) / (vmax * (a$ar/sr) - (a$af * a$nf + a$aw * a$nw + a$ar * a$nr))
    
    NPP_NC <- (U0 - nleach) / (nwood + nburial)   # will be in g C m-2 yr-1
    NPP_N <- NPP_NC*10^-3 # returned in kg C m-2 yr-1
    
    df <- data.frame(NPP_N, nwood,nburial,nleach,a$aw)
    return(df)   
}