#### To compute autotrophic respiration
Compute_Ra <- function(a, NPP) {
    
    # parameters
    # a is variable
    # a1 = intercept for Rleaf and Leaf N;
    # b1 = exponent for Rleaf and Leaf N;
    # a2 = intercept for Rstem and Stem N;
    # b2 = exponent for Rstem and Stem N;
    # a3 = intercept for Rroot and Root N;
    # b3 = exponent for Rroot and Root N;
    # rleaf: leaf respiration, kg m-2 yr-1
    # rstem: stem respiration, kg m-2 yr-1
    # rroot: root respiration, kg m-2 yr-1
    # conv: unit conversion for MOL_N_TO_GRAMS_N * MOL_2_MMOL;
    # sf: leaf turnover time yr-1
    # ss: stem turnover time yr-1
    # sr: root turnover time yr-1
    
    ## Set unit conversion factor;
    MOL_N_TO_GRAMS_N  <- 14.0
    MOL_2_MMOL <- 1000.0
    
    ## Use shoot and root nc ratio to obtain mmol [N] g-1
    leafn <- a$nf / MOL_N_TO_GRAMS_N * MOL_2_MMOL
    stemn <- a$nw / MOL_N_TO_GRAMS_N * MOL_2_MMOL
    rootn <- a$nr / MOL_N_TO_GRAMS_N * MOL_2_MMOL
    
    ## calculate dark respiration for leaf, stem and root, (nmol g-1 s-1)
    respl <- a1 * (leafn^b1)
    resps <- a2 * (stemn^b2)
    respr <- a3 * (rootn^b3)
    
    ## unit conversion - per kg C allocated to leaf, stem and root (nmol m-2 s-1)
    respl <- respl * (NPP * 1000.0 * a$af / sf / cfrac) # need to multiple by biomass, not npp
    resps <- resps * (NPP * 1000.0 * a$aw / ss / cfrac)
    respr <- respr * (NPP * 1000.0 * a$ar / sr / cfrac)
    
    ## unit conversion, to kg C m-2 yr-1
    rleaf <- respl * 3600.0 * 24.0 * 365.0 * 10^-9 * 0.012
    rstem <- resps * 3600.0 * 24.0 * 365.0 * 10^-9 * 0.012
    rroot <- respr * 3600.0 * 24.0 * 365.0 * 10^-9 * 0.012
    
    ## Compute autotrophic respiration
    r_autotrophic <- rleaf + rstem + rroot
    
    return(r_autotrophic)
    
    
}

#### To compute autotrophic respiration
Compute_Rdark <- function(nfdf, pfdf, NPP) {
    
    # parameters
    # a is variable
    # a1 = intercept for Rleaf and Leaf N;
    # b1 = exponent for Rleaf and Leaf N;
    # a2 = intercept for Rstem and Stem N;
    # b2 = exponent for Rstem and Stem N;
    # a3 = intercept for Rroot and Root N;
    # b3 = exponent for Rroot and Root N;
    # rleaf: leaf respiration, kg m-2 yr-1
    # rstem: stem respiration, kg m-2 yr-1
    # rroot: root respiration, kg m-2 yr-1
    # conv: unit conversion for MOL_N_TO_GRAMS_N * MOL_2_MMOL;
    # sf: leaf turnover time yr-1
    # ss: stem turnover time yr-1
    # sr: root turnover time yr-1
    
    ## Set unit conversion factor;
    MOL_N_TO_GRAMS_N  <- 14.0
    MOL_2_MMOL <- 1000.0
    
    ## Temperature of warmest quarter
    twq <- tk - 273.5
    
    ## calculate jmax and vcmax
    ncontent <- nfdf$nf * NPP * nfdf$af / sf    # g N m-2
    pcontent <- pfdf$pf * NPP * pfdf$af / sf
    
    # convert SLA from m2 kg-1 to m2 g-1
    sla_m2_per_g <- SLA / 1000.0
    
    N0 <- ncontent * kn / (1.0 - exp(-kn * sla_m2_per_g*nfdf$af*NPP/sf/cfrac))
    P0 <- pcontent * kn / (1.0 - exp(-kn * sla_m2_per_g*pfdf$af*NPP/sf/cfrac))
    
    log_vcmax <- 3.946 + 0.921 * log(N0) + 0.121 * log(P0) + 0.282 * log(N0) * log(P0)
    vcmax <- exp(log_vcmax)

    vcmax25 <- vcmax
    
    #browser()
    
    ## calculate rdark leaf at 25 C
    r_leaf_dark_25 <- 1.2636 + 0.0728 * N0 + 0.015 * P0 + 0.0095 * vcmax25 - 0.0358 * twq
    
    # convert unit from umol CO2 m-2 s-1 to kg C m-2 yr-1
    r_leaf_dark_25_out <- r_leaf_dark_25 * 1E-6 * 12.0 * 3600 * 24 * 365 / 1000.0
    
    return(r_leaf_dark_25_out)
}



