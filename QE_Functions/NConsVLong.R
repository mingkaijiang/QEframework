### Calculate the very long term nutrient cycling constraint for N, i.e. passive pool equilibrated
# it is just Nin = Nleach
NConsVLong <- function(df, a) {
    # passed are bf and nf, the allocation and plant N:C ratios
    # parameters : 
    # Nin is fixed N inputs (N fixed and deposition) in g m-2 yr-1 (could vary fixation)
    # leachn is the rate of leaching of the mineral N pool (per year)
    
    # equation for N constraint with just leaching
    U0 <- Nin
    nleach <- leachn/(1-leachn) * (a$nfl*a$af + a$nr*(a$ar) + a$nw*a$aw)
    NPP_NC <- U0 / (nleach)   # will be in g C m-2 yr-1
    NPP_N <- NPP_NC*10^-3     # returned in kg C m-2 yr-1
    
    df <- data.frame(NPP_N,nleach)
    return(df)   
}

### Calculate the very long term nutrient cycling constraint for N, i.e. passive pool equilibrated
# it is just Nin = Nleach
# specifically for explicit mineral N pool
NConsVLong_expl_min <- function(df, a) {
    # passed are bf and nf, the allocation and plant N:C ratios
    # parameters : 
    # Nin is fixed N inputs (N fixed and deposition) in g m-2 yr-1 (could vary fixation)
    # leachn is the rate of leaching of the mineral N pool (per year)
    # nuptakerate is the rate of N uptake [yr-1] 
    # Nmin is the mineral N pool
    
    # equation for N constraint with just leaching
    U0 <- Nin
    nleach <- leachn
    
    Nmin <- U0 / nleach
    NPP_NC <- Nmin * nuptakerate / (a$nfl*a$af + a$nr*a$ar + a$nw*a$aw)
    # The above equation is the same as:
    # NPP_NC <- U0 * nuptakerate / (nleach * (a$nfl*a$af + a$nr*(a$ar) + a$nw*a$aw))
    NPP_N <- NPP_NC*10^-3     # returned in kg C m-2 yr-1
    
    df <- data.frame(NPP_N,nleach)
    return(df)   
}

### Calculate the very long term nutrient cycling constraint for N, i.e. passive pool equilibrated
# it is just Nin = Nleach
# specifically for nuptake as a function of root biomass - O-CN approach
# i.e. N uptake as a saturating function of mineral N
NConsVLong_root_ocn <- function(carb_diox) {                   
    # passed are bf and nf, the allocation and plant N:C ratios
    # parameters : 
    # Nin is fixed N inputs (N fixed and deposition) in g m-2 yr-1 (could vary fixation)
    # leachn is the rate of leaching of the mineral N pool (per year)
    # nuptakerate is the rate of N uptake [yr-1] 
    # Nmin is the mineral N pool
    # sr is the decay rate of root in yr-1
    # k - empirically dervied for now
    # vmax
    
    # allocation coefficients
    ar <- aroot
    af <- aleaf
    aw <- 1 - ar - af
    
    # compute Nmin
    Nmin <- Nin / leachn
    
    # arguments
    arg1 <- (ar / sr) * vmax
    arg2 <- Nmin / (Nmin + k)
    arg3 <- arg1 * arg2
    
    # N concentrations of rest of plant   # in g N g-1 C
    if (nwvar == FALSE) {
        nw <- nwood
        nf <- (arg3 - nw * aw) / (af + nrho * ar)
    } else {
        nf <- arg3 / (af + nrho * ar + nwood * aw)
    }
    
    # allocation coefficients
    a_nf <- allocn(nf)
    
    # solve for equilibrium npp
    npp <- photo_constraint_full_cn(nf,a_nf,carb_diox)
    
    ans <- data.frame(nf,npp)
    colnames(ans) <- c("equilnf","equilNPP")
    
    return(ans)   
}

NConsVLong_root_ocn_original <- function(df, a) {                   
    # passed are bf and nf, the allocation and plant N:C ratios
    # parameters : 
    # Nin is fixed N inputs (N fixed and deposition) in g m-2 yr-1 (could vary fixation)
    # leachn is the rate of leaching of the mineral N pool (per year)
    # nuptakerate is the rate of N uptake [yr-1] 
    # Nmin is the mineral N pool
    # sr is the decay rate of root in yr-1
    # k - empirically dervied for now
    # vmax
    
    # compute N mineral pool
    Nmin <- k * (a$nfl*a$af + a$nr*a$ar + a$nw*a$aw) / (a$ar / sr - (a$nfl*a$af + a$nr*a$ar + a$nw*a$aw))
    
    browser()
    
    # equation for N constraint with just leaching
    U0 <- Nin
    nleach <- leachn * Nmin
    
    NPP_NC <- U0 / nleach
    NPP_N <- NPP_NC*10^-3     # returned in kg C m-2 yr-1
    
    df <- data.frame(NPP_N,nleach)
    return(df)   
}

### Calculate the very long term nutrient cycling constraint for N, i.e. passive pool equilibrated
# it is just Nin = Nleach
# specifically for nuptake as a function of root biomass - GDAY approach
# i.e. N uptake as a saturating function of root biomass
NConsVLong_root_gday <- function(df, a) {                   # why this small?
    # passed are bf and nf, the allocation and plant N:C ratios
    # parameters : 
    # Nin is fixed N inputs (N fixed and deposition) in g m-2 yr-1 (could vary fixation)
    # leachn is the rate of leaching of the mineral N pool (per year)
    # nuptakerate is the rate of N uptake [yr-1] 
    # Nmin is the mineral N pool
    # sr is the decay rate of root in yr-1
    # kr is the value of root carbon at which 50% of the available N is taken up  
    
    
    # equation for N constraint with just leaching
    U0 <- Nin
    Nmin <- Nin / leachn
    A_NF <- a$nfl*a$af + a$nr*a$ar + a$nw*a$aw
    root_biomass <- a$ar / sr
    nleach <- Nmin * leachn
    
    #browser()
    
    # equation for NPP
    #NPP_NC <- (root_biomass * Nmin - (A_NF * kr)) / (A_NF * root_biomass)
    NPP_NC <- (Nmin - kr * A_NF / root_biomass ) / A_NF
    NPP_N <- NPP_NC*10^-3     # returned in kg C m-2 yr-1
    
    #browser()
    
    df <- data.frame(NPP_N, nleach)
    return(df)   
}


