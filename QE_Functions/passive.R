### Burial fractions from passive pool
passive <- function(df, a) {
    
    len <- length(df)
    
    decomp <- 0.00013*52*Actsoil(Tsoil)   # decomposition of passive pool per year without priming
    
    # re-burial fraction = fraction of C released from passive pool that is re-buried in it
    pas <- 0.996 - (0.85-0.68*Texture)
    psa <- 0.42
    ppa <- 0.45
    pap <- 0.004
    psp <- 0.03
    qq <-  ppa*(pap + psp*pas)/(1-pas*psa)   # re-burial fraction
    
    muf <- c()
    mur <- c()
    omegaf <- c()
    omegar <- c()
    transfer_fa <- c()
    transfer_ra <- c()
    
    
    # transfer coefficients among litter and soil pools
    for (i in 1:len) {
        muf[i] <- max(0,min(0.85 - 0.018*ligfl/cfrac/a[i, "nfl"],1))    
        mur[i] <- max(0,min(0.85 - 0.018*ligrl/cfrac/a[i, "nr"],1))
    }
    pma <- pna <- 0.45
    pua <- 0.55*(1-ligfl)
    pus <- 0.7*ligfl
    pva <- 0.45*(1-ligrl)
    pvs <- 0.7*ligrl
    
    # burial fractions for foliage (omegaf) and root (omegar) into passive pool
    det <- 1-psa*pas
    omegau <- (pap*(pua+pus*psa+psp*(pus+pua*pas)))/det
    omegav <- (pap*(pva+pus*psa+psp*(pvs+pva*pas)))/det
    omegam <- (pap*pma+psp*pma*pas)/det
    omegaf <- muf*omegam + (1-muf)*omegau
    omegar <- mur*omegam + (1-mur)*omegav  
    
    # fraction of foliage and root litter being transferred to active pool
    transfer_fa <- muf*pma + (1-muf)*psa
    transfer_ra <- mur*pma + (1-mur)*psa
    
    ret <- data.frame(decomp, qq, omegaf, omegar, transfer_fa, transfer_ra)
    
    return(ret)
}

### Burial fractions from passive pool
# specifically for exudation
passive_exudation <- function(df, a, npp) {
    
    len <- length(df)
    
    # re-burial fraction = fraction of C released from passive pool that is re-buried in it
    pas <- 0.996 - (0.85-0.68*Texture)       # e.g. from active to slow
    psa <- 0.42     # 0.42; 0.396 (slow)
    ppa <- 0.45     # should be a function of exudate C
    pap <- 0.004
    psp <- 0.03     # 0.03; 0.028 (slow)
    qq <-  ppa*(pap + psp*pas)/(1-pas*psa)   # re-burial fraction of the passive pool, consider to adjust
    
    #browser()
    
    muf <- c()
    mur <- c()
    omegaf <- c()
    omegar <- c()
    transfer_fa <- c()
    transfer_ra <- c()
    
    # transfer coefficients among litter and soil pools
    for (i in 1:len) {
        muf[i] <- max(0,min(0.85 - 0.018*ligfl/cfrac/a[i, "nfl"],1))    
        mur[i] <- max(0,min(0.85 - 0.018*ligrl/cfrac/a[i, "nr"],1))
    }
    pma <- pna <- 0.45
    pua <- 0.55*(1-ligfl)
    pus <- 0.7*ligfl
    pva <- 0.45*(1-ligrl)
    pvs <- 0.7*ligrl
    
    # burial fractions for foliage (omegaf) and root (omegar) into passive pool
    det <- 1-psa*pas
    omegau <- (pap*(pua+pus*psa+psp*(pus+pua*pas)))/det
    omegav <- (pap*(pva+pus*psa+psp*(pvs+pva*pas)))/det
    omegam <- (pap*pma+psp*pma*pas)/det
    omegaf <- muf*omegam + (1-muf)*omegau
    omegar <- mur*omegam + (1-mur)*omegav  
    
    # fraction of foliage and root litter being transferred to active pool
    transfer_fa <- muf*pma + (1-muf)*psa
    transfer_ra <- mur*pma + (1-mur)*psa
    
    # total active out
    active_plant_in <- transfer_fa * (npp * a$af / sf) + transfer_ra * (npp + (a$ar-a$ar*a$ariz) / sr)
    c_into_exud <- npp * a$ar * a$ariz
    active_tot_in <- active_plant_in + c_into_exud 
    
    # decomposition rate of the passive pool
    decomp <- adjust_passive_residence_time(df, a, active_tot_in) * Actsoil(Tsoil)   

    ret <- data.frame(decomp, qq, omegaf, omegar, transfer_fa, transfer_ra)
    
    return(ret)
}
