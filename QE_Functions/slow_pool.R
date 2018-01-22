### Burial fractions from passive pool
slow_pool <- function(df, a) {
    
    len <- length(df)
    
    decomp_p <- 0.00013*52*Actsoil(Tsoil)   # decomposition of passive pool per year without priming
    decomp_s <- 0.0038*52*Actsoil(Tsoil)    # decomposition of slow pool per year without priming
    
    # re-burial fraction 
    pas <- 0.996 - (0.85-0.68*Texture)
    psa <- 0.42
    ppa <- 0.45
    pap <- 0.004
    psp <- 0.03
    
    # fraction of C released from passive pool that is re-buried in it
    qpq <-  ppa*(pap + psp*pas)/(1-pas*psa)   
    # fraction of C released from slow pool that is re-buried in it
    qsq <- (pas*ppa*psp + pas*psa) / (1 - ppa*pap)
    
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
    omegafp <- muf*omegam + (1-muf)*omegau
    omegarp <- mur*omegam + (1-mur)*omegav  
    
    # burial fractions for foliage and root into slow pool
    omegaus <- (pus+pas*pua)/det
    omegavs <- (pvs+pas*pva)/det
    omegams <- (pas*pma)/det
    omegafs <- muf*omegams + (1-muf)*omegaus
    omegars <- mur*omegams + (1-mur)*omegavs
    
    # fraction of foliage and root litter being transferred to active pool
    transfer_fa <- muf*pma + (1-muf)*psa
    transfer_ra <- mur*pma + (1-mur)*psa
    
    ret <- data.frame(decomp_p, decomp_s, qpq, qsq,
                      omegafp, omegarp, omegafs, omegars, 
                      transfer_fa, transfer_ra)
    
    return(ret)
}
