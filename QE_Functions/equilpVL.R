
# Find very-long term equilibrated pf based on equilibrated NPP calculated from equilnf profile
equilpVL <- function(equilNPP) {

    ar <- aroot
    af <- aleaf
    aw <- 1 - ar - af
    
    # prepare very long term phosphorus fluxes
    U0 = Pin
    pleach <- leachp/(1-leachp-k1)
    pocc <- (k3/(k2+k3))*(k1/(1-k1-pleach))
    
    # Convert NPP unit
    NPP_PC <- equilNPP*10^3     # convert to g C m-2 yr-1
    
    # Calculate equilnf, based on equilNPP_P
    Y <- U0/(NPP_PC*(pleach+pocc))
    if(pwvar == FALSE) {
        pf <- (Y - pwood * aw) / ((1.0 - pretrans) * af + prho * ar)
    } else {
        pf <- Y / ((1.0 - pretrans) * af + prho * ar + aw * pwood)
    }
    
    # obtain equilpf  
    return(pf)
}


# Find very-long term equilibrated pf based on equilibrated NPP calculated from equilnf profile
# specifically for explicit mineral N and P pool
equilpVL_expl_min <- function(equilNPP, Pin = 0.02, leachp=0.05,
                              pwvar = TRUE, pwood = 0.0003, prho = 0.7,
                              pretrans = 0.6, k1 = 0.01, k2 = 0.01, k3 = 0.05, puptakerate = 0.82395) {
    # prepare allocation partitioning
    ar <- 0.2
    af <- 0.2
    aw <- 1 - ar - af
    
    # prepare very long term phosphorus fluxes
    U0 = Pin
    pleach <- leachp/(1-leachp-k1)
    pocc <- (k3/(k2+k3))*(k1/(1-k1-pleach))
    
    # Convert NPP unit
    NPP_PC <- equilNPP*10^3     # convert to g C m-2 yr-1
    
    # Calculate equilnf, based on equilNPP_P
    Y <- U0 / (NPP_PC*(pleach+pocc))
    if(pwvar == FALSE) {
        pf <- (Y - pwood * aw) / ((1.0 - pretrans) * af + prho * ar)
    } else {
        pf <- Y / ((1.0 - pretrans) * af + prho * ar + aw * pwood)
    }
    
    # obtain equilpf  
    return(pf)
}

# Find very-long term equilibrated pf based on equilibrated NPP calculated from equilnf profile
# specifically for N uptake as a function of biomass - OCN approach
# i.e. N uptake as a saturating function of mineral N
equilpVL_root_ocn <- function(equilNPP, Pin = 0.02, leachp=0.05,
                              pwvar =F, pwood = 0.0003, prho = 0.7,
                              pretrans = 0.6, k1 = 0.01, k2 = 0.01, k3 = 0.05, 
                              puptakerate = 0.82395) {
    # prepare allocation partitioning
    ar <- 0.2
    af <- 0.2
    aw <- 1 - ar - af
    
    # prepare very long term phosphorus fluxes
    U0 = Pin
    pleach <- leachp/(1-leachp-k1)
    pocc <- (k3/(k2+k3))*(k1/(1-k1-pleach))
    
    # Convert NPP unit
    NPP_PC <- equilNPP*10^3     # convert to g C m-2 yr-1
    
    # Calculate equilnf, based on equilNPP_P
    Y <- U0 / (NPP_PC*(pleach+pocc))
    if(pwvar == FALSE) {
        pf <- (Y - pwood * aw) / ((1.0 - pretrans) * af + prho * ar)
    } else {
        pf <- Y / ((1.0 - pretrans) * af + prho * ar + aw * pwood)
    }
    
    # obtain equilpf  
    return(pf)
}

# Find very-long term equilibrated pf based on equilibrated NPP calculated from equilnf profile
# specifically for N uptake as a function of biomass - GDAY approach
# i.e. N uptake as a saturating function of root biomass
equilpVL_root_gday <- function(equilNPP, Pin = 0.02, leachp=0.05,
                              pwvar = TRUE, pwood = 0.0003, prho = 0.7,
                              pretrans = 0.6, k1 = 0.01, k2 = 0.01, k3 = 0.05, 
                              puptakerate = 0.82395) {
    # prepare allocation partitioning
    ar <- 0.2
    af <- 0.2
    aw <- 1 - ar - af
    
    # prepare very long term phosphorus fluxes
    U0 = Pin
    pleach <- leachp/(1-leachp-k1)
    pocc <- (k3/(k2+k3))*(k1/(1-k1-pleach))
    
    # Convert NPP unit
    NPP_PC <- equilNPP*10^3     # convert to g C m-2 yr-1
    
    # Calculate equilnf, based on equilNPP_P
    Y <- U0 / (NPP_PC*(pleach+pocc))
    if(pwvar == FALSE) {
        pf <- (Y - pwood * aw) / ((1.0 - pretrans) * af + prho * ar)
    } else {
        pf <- Y / ((1.0 - pretrans) * af + prho * ar + aw * pwood)
    }
    
    # obtain equilpf  
    return(pf)
}


# Find very-long term equilibrated pf based on equilibrated NPP calculated from equilnf profile
# specifically for exudation
equilpVL_exudation <- function(equilNPP, Pin = 0.02, leachp=0.05,
                              pwvar = TRUE, pwood = 0.0003, prho = 0.7,
                              pretrans = 0.6, k1 = 0.01, k2 = 0.01, k3 = 0.05) {
    # prepare allocation partitioning
    ar <- 0.15
    af <- 0.2
    ae <- 0.05
    aw <- 1 - ar - af - ae
    
    # prepare very long term phosphorus fluxes
    U0 = Pin
    pleach <- leachp/(1-leachp-k1)
    pocc <- (k3/(k2+k3))*(k1/(1-k1-pleach))
    
    # Convert NPP unit
    NPP_PC <- equilNPP*10^3     # convert to g C m-2 yr-1
    
    # Calculate equilnf, based on equilNPP_P
    Y <- U0/(NPP_PC*(pleach+pocc))
    if(pwvar == FALSE) {
        pf <- (Y - pwood * aw) / ((1.0 - pretrans) * af + prho * ar)
    } else {
        pf <- Y / ((1.0 - pretrans) * af + prho * ar + aw * pwood)
    }
    
    # obtain equilpf  
    return(pf)
}