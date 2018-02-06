


### Function for nutrient N constraint in longterm ie passive, leaching, wood and CWD considered
NConsLong_CWD <- function(df, a, Nin=0.4, leachn=0.05, 
                          Tsoil = 15, Texture = 0.5, ligfl = 0.2, ligrl = 0.16,
                          Cpass = 2680, ncp = 0.1, sw = 0.02) {
    # passed are df and a, the allocation and plant N:C ratios
    # parameters : 
    # Nin is fixed N inputs (N deposition annd fixation) in g m-2 yr-1 (could vary fixation)
    # nleach is the rate of n leaching of the mineral pool (per year)
    # Tsoil is effective soil temperature for decomposition
    # Texture is the fine soil fraction
    # ligfl and ligrl are the lignin:C fractions in the foliage and root litter
    # Cpass is the passive pool size in g C m-2
    # ncp is the NC ratio of the passive pool in g N g-1 C
    # sw is the decay rate of wood in yr-1
    # ncwd is the rate of N storage in CWD pool in g N m-2 yr-1
    
    # passive pool burial 
    pass <- passive(df, a, Tsoil, Texture, ligfl, ligrl)
    omegap <- a$af*pass$omegaf + a$ar*pass$omegar 
    
    # equation for N constraint with passive, wood, CWD and leaching
    U0 <- Nin + (1-pass$qq) * pass$decomp * Cpass * ncp   # will be a constant if decomp rate is constant
    nwood <- a$aw*a$nw
    nburial <- omegap*ncp
    nleach <- leachn/(1-leachn) * (a$nfl*a$af + a$nr*a$ar + a$nw*a$aw)
    ncwd <- nwood*sw   # do we need to consider fluxes lost from this pool? e.g. into active and slow?

    NPP_NC <- U0 / (ncwd + nburial + nleach)   # will be in g C m-2 yr-1
    NPP_N <- NPP_NC*10^-3 # returned in kg C m-2 yr-1
    
    df <- data.frame(NPP_N, nwood,nburial,nleach,a$aw)
    return(df)   
}


### Function for nutrient N constraint in longterm ie passive, leaching, wood considered
# specifically for explicit mineral N
NConsLong_expl_min <- function(df, a, Cpass, NinL) {
    # passed are df and a, the allocation and plant N:C ratios
    # parameters : 
    # Nin is fixed N inputs (N deposition annd fixation) in g m-2 yr-1 (could vary fixation)
    # nleach is the rate of n leaching of the mineral pool (per year)
    # Tsoil is effective soil temperature for decomposition
    # Texture is the fine soil fraction
    # ligfl and ligrl are the lignin:C fractions in the foliage and root litter
    # Cpass is the passive pool size in g C m-2
    # ncp is the NC ratio of the passive pool in g N g-1 C
    # nuptakerate is the rate of N uptake in yr-1
    
    # passive pool burial 
    pass <- passive(df, a)
    omegap <- a$af*pass$omegaf + a$ar*pass$omegar 
    
    # equation for N constraint with passive, wood, and leaching
    U0 <- NinL + (1-pass$qq) * pass$decomp * Cpass * ncp   # will be a constant if decomp rate is constant
    nwood <- 0 # a$aw*a$nw
    nburial <- omegap*ncp

    #NPP_NC <- ((U0 - nwood - nburial) / leachn) * nuptakerate / (a$nfl*a$af + a$nr*a$ar + a$nw*a$aw)   # will be in g C m-2 yr-1
    NPP_NC <- (nuptakerate * U0) / ((a$nfl*a$af + a$nr*a$ar + a$nw*a$aw) * leachn + nuptakerate * nwood + nuptakerate * nburial)
    
    NPP <- NPP_NC*10^-3 # returned in kg C m-2 yr-1
    
    nleach <- leachn * (NPP_NC * (a$nfl*a$af + a$nr*a$ar + a$nw*a$aw)) /nuptakerate
    
    df <- data.frame(NPP, nwood,nburial,nleach,a$aw)
    return(df)   
}

### Function for nutrient N constraint in longterm ie passive, leaching, wood considered
# specifically for uptake as a function of root  - O-CN approach
# i.e. N uptake as a saturating function of mineral N
NConsLong_root_ocn <- function(df, a, Cpass, NinL) {
    # passed are df and a, the allocation and plant N:C ratios
    # parameters : 
    # Nin is fixed N inputs (N deposition annd fixation) in g m-2 yr-1 (could vary fixation)
    # nleach is the rate of n leaching of the mineral pool (per year)
    # Tsoil is effective soil temperature for decomposition
    # Texture is the fine soil fraction
    # ligfl and ligrl are the lignin:C fractions in the foliage and root litter
    # Cpass is the passive pool size in g C m-2
    # ncp is the NC ratio of the passive pool in g N g-1 C
    # nuptakerate is the rate of N uptake in yr-1
    # sr is the decay rate of root in yr-1
    # k - empirically derived
    # vmax in real O-CN vmax down, leaf NC up
    # NPP = (Nin - λloss A K nf / (ar / sr Vmax – A nf)) / Ωp (nf) np
    
    # passive pool burial 
    pass <- passive(df, a)
    omegap <- a$af*pass$omegaf + a$ar*pass$omegar 
    
    # equation for N constraint with passive, wood, and leaching
    U0 <- NinL + (1-pass$qq) * pass$decomp * Cpass * ncp   # will be a constant if decomp rate is constant
    nwood <- 0 # a$aw*a$nw
    nburial <- omegap*ncp
    nleach <- leachn * k * (a$af * a$nf + a$aw * a$nw + a$ar * a$nr) / (vmax * (a$ar/sr) - (a$af * a$nf + a$aw * a$nw + a$ar * a$nr))
    
    NPP_NC <- (U0 - nleach) / (nwood + nburial)   # will be in g C m-2 yr-1
    NPP_N <- NPP_NC*10^-3 # returned in kg C m-2 yr-1
    
    #browser()
    
    df <- data.frame(NPP_N, nwood,nburial,nleach,a$aw)
    return(df)   
}

### Function for nutrient N constraint in longterm ie passive, leaching, wood considered
# specifically for uptake as a function of root - GDAY approach
# i.e. N uptake as a saturating function of root biomass
NConsLong_root_gday <- function(df, a, NinL, 
                               Cpass) {
    # passed are df and a, the allocation and plant N:C ratios
    # parameters : 
    # Nin is fixed N inputs (N deposition annd fixation) in g m-2 yr-1 (could vary fixation)
    # nleach is the rate of n leaching of the mineral pool (per year)
    # Tsoil is effective soil temperature for decomposition
    # Texture is the fine soil fraction
    # ligfl and ligrl are the lignin:C fractions in the foliage and root litter
    # Cpass is the passive pool size in g C m-2
    # ncp is the NC ratio of the passive pool in g N g-1 C
    # nuptakerate is the rate of N uptake in yr-1
    # sr is the decay rate of root in yr-1
    # kr is the value of root carbon at which 50% of the available N is taken up
    
    # passive pool burial 
    pass <- passive(df, a)
    omegap <- a$af*pass$omegaf + a$ar*pass$omegar 
    
    # equation for N constraint with passive, wood, and leaching
    U0 <- NinL + (1-pass$qq) * pass$decomp * Cpass * ncp   # will be a constant if decomp rate is constant
    nwood <- 0 # a$aw*a$nw
    nburial <- omegap*ncp
    
    A_NF <- a$nfl*a$af + a$nr*a$ar + a$nw*a$aw

    # equation for NPP
    NPP_NC <- (U0 - A_NF * leachn * kr * (sr / a$ar)) / (nwood + nburial + A_NF*leachn)
    
    NPP_N <- NPP_NC*10^-3 # returned in kg C m-2 yr-1
    
    nleach <- leachn * (U0 - nwood * NPP_N - nburial * NPP_N)
    
    df <- data.frame(NPP_N, nwood,nburial,nleach,a$aw)
    return(df)   
}



### Function for nutrient N constraint in longterm ie passive, leaching, wood considered
### specifically for the case of passive NC ratio depends on soil mineral N
NConsLong_variable_pass <- function(df, a, Nin=0.4, leachn=0.05, eqNPP,
                                    Tsoil = 15, Texture = 0.5, ligfl = 0.2, ligrl = 0.16,
                                    Cpass = 2680, nuptakerate = 0.96884) {
    # passed are df and a, the allocation and plant N:C ratios
    # parameters : 
    # Nin is fixed N inputs (N deposition annd fixation) in g m-2 yr-1 (could vary fixation)
    # nleach is the rate of n leaching of the mineral pool (per year)
    # Tsoil is effective soil temperature for decomposition
    # Texture is the fine soil fraction
    # ligfl and ligrl are the lignin:C fractions in the foliage and root litter
    # Cpass is the passive pool size in g C m-2
    # ncp is the NC ratio of the passive pool in g N g-1 C
    # nup is the rate of nuptake in /yr
    # n1 is the intercept of function U = n1 + n2 * Nmin
    # n2 is the coefficient of function U = n1 + n2 * Nmin

    # passive pool burial 
    pass <- passive(df, a, Tsoil, Texture, ligfl, ligrl)
    omegap <- a$af*pass$omegaf + a$ar*pass$omegar 
    
    # calculate ncp (a linear function of mineral N)
    ncp <- calc_passive_nc(equilNPP=eqNPP, adf = a, nuptake = nuptakerate)
    
    # equation for N constraint with passive, wood, and leaching
    U0 <- Nin + (1-pass$qq) * pass$decomp * Cpass * ncp   # will be a constant if decomp rate is constant
    nwood <- a$aw*a$nw
    nburial <- omegap*ncp
    nleach <- leachn/(1-leachn) * (a$nfl*a$af + a$nr*(a$ar) + a$nw*a$aw)
    
    NPP_NC <- U0 * nuptakerate / (nwood + nburial + nleach)   # will be in g C m-2 yr-1
    
    NPP_N <- NPP_NC*10^-3 # returned in kg C m-2 yr-1
    
    df <- data.frame(NPP_N, nwood,nburial,nleach,a$aw)
    return(df)   
}


### Function for nutrient N constraint in longterm ie passive, leaching, wood considered
# specifically for exudation
NConsLong_exudation <- function(df, a, Cpass, NinL) {
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
    
    for (i in 1:len) {
        fPC <- function(NPP) {
            # passive pool burial 
            pass <- passive_exudation(df[i], a[i,], NPP)   
            omegap <- a[i,]$af*pass$omegaf + a[i,]$ar*pass$omegar 
            
            # equation for N constraint with passive, wood, and leaching
            U0 <- NinL + (1-pass$qq) * pass$decomp * Cpass * ncp   
            nwood <- 0 # a[i,]$aw*a[i,]$nw
            nburial <- omegap*ncp
            nleach <- leachn / (1-leachn) * (a[i,]$nfl*a[i,]$af + a[i,]$nr*a[i,]$ar 
                                             + a[i,]$nw*a[i,]$aw)
            
            NPP_NC <- U0 / (nleach + nwood + nburial)
            
            out <- NPP_NC*10^-3 - NPP # returned in kg C m-2 yr-1
        }
        
        ans[i] <- uniroot(fPC,interval=c(0.1,20), trace=T)$root
        
    }

    out <- data.frame(ans, a$aw)
    colnames(out) <- c("NPP", "aw")
    return(out)   
}

### Function for nutrient N constraint in longterm ie passive, leaching, wood considered
# specifically for exudation
NConsLong_exudation_medium <- function(df, a, Cpass, NinL) {
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
    
    for (i in 1:len) {
        fPC <- function(NPP) {
            # passive pool burial 
            pass <- slow_exudation(df[i], a[i,], NPP)   
            omegap <- a[i,]$af*pass$omegafp + a[i,]$ar*pass$omegarp 
            
            # equation for N constraint with passive, wood, and leaching
            U0 <- NinL + (1-pass$qpq) * pass$decomp_p * Cpass * ncp   
            nwood <- 0 # a[i,]$aw*a[i,]$nw
            nburial <- omegap*ncp
            nleach <- leachn / (1-leachn) * (a[i,]$nfl*a[i,]$af + a[i,]$nr*a[i,]$ar 
                                             + a[i,]$nw*a[i,]$aw)
            
            NPP_NC <- U0 / (nleach + nwood + nburial)
            
            out <- NPP_NC*10^-3 - NPP # returned in kg C m-2 yr-1
        }
        
        ans[i] <- uniroot(fPC,interval=c(0.1,20), trace=T)$root
        
    }
    
    out <- data.frame(ans, a$aw)
    colnames(out) <- c("NPP", "aw")
    return(out)   
}
