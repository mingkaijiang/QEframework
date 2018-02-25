### Burial fractions from passive pool
passive <- function(nf, a, Tsoil = 15, Texture = 0.5, ligfl = 0.2, ligrl = 0.16) {
    
    decomp <- 0.00013*52*Actsoil(Tsoil)   # decomposition of passive pool per year without priming
    
    # re-burial fraction = fraction of C released from passive pool that is re-buried in it
    pas <- 0.996 - (0.85-0.68*Texture)
    psa <- 0.42
    ppa <- 0.45
    pap <- 0.004
    psp <- 0.03
    qq <-  ppa*(pap + psp*pas)/(1-pas*psa)   # re-burial fraction
    
    # transfer coefficients among litter and soil pools
    cfrac <- 0.45
    
    muf <- c()
    mur <- c()
    
    for (i in 1:length(a$nfl)) {
        muf[i] <- max(0,min(0.85 - 0.018*ligfl/cfrac/a$nfl[i],1))   
        mur[i] <- max(0,min(0.85 - 0.018*ligrl/cfrac/a$nr[i],1))    
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

