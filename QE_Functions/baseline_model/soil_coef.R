### Burial fractions from passive pool
soil_coef <- function(df, a) {
    # Texture: soil texture value - pre-defined
    # Tsoil: soil temperature
    
    ### get the dataframe length
    len <- length(df)
    
    ### calculate decomposition parameters
    ## decomposition of passive pool per year
    decomp_pass <- 0.00013*52*soil_decomp(Tsoil) 
    
    ## decomposition of slow pool per year
    ## decomp of slow pool is 30 times faster than that of the passive pool
    ## is it too fast? 15 year turnover!
    decomp_slow <- 0.0038*52*soil_decomp(Tsoil)    
    
    ### re-burial fraction 
    pas <- 0.996 - (0.85-0.68*Texture)
    psa <- 0.42
    ppa <- 0.45
    pap <- 0.004
    psp <- 0.03
    
    ### fraction of C released from passive pool that is re-buried in it
    qq_pass <-  ppa*(pap + psp*pas)/(1-pas*psa)
    
    ### fraction of C released from slow pool that is re-buried in it
    qq_slow <- (pas*ppa*psp + pas*psa) / (1 - ppa*pap)
    
    ### prepare storage df
    mu_f <- c()
    mu_r <- c()
    omega_af_pass <- c()
    omega_ar_pass <- c()
    omega_af_slow <- c()
    omega_ar_slow <- c()
    transfer_fa <- c()
    transfer_ra <- c()
    
    
    ### transfer coefficients among litter and soil pools
    for (i in 1:len) {
        mu_f[i] <- max(0,min(0.85 - 0.018*ligfl/cfrac/a[i, "nfl"],1))    
        mu_r[i] <- max(0,min(0.85 - 0.018*ligrl/cfrac/a[i, "nr"],1))
    }
    pma <- pna <- 0.45
    pua <- 0.55*(1-ligfl)
    pus <- 0.7*ligfl
    pva <- 0.45*(1-ligrl)
    pvs <- 0.7*ligrl
    
    ### burial fractions for foliage (omegaf) and root (omegar) into passive pool
    det <- 1-psa*pas
    omega_au_pass <- (pap*(pua+pus*psa+psp*(pus+pua*pas)))/det
    omega_av_pass <- (pap*(pva+pus*psa+psp*(pvs+pva*pas)))/det
    omega_am_pass <- (pap*pma+psp*pma*pas)/det
    omega_af_pass <- mu_f*omega_am_pass + (1-mu_f)*omega_au_pass
    omega_ar_pass <- mu_r*omega_am_pass + (1-mu_r)*omega_av_pass  
    
    ### burial fractions for foliage and root into slow pool
    omega_au_slow <- (pus+pas*pua)/det
    omega_av_slow <- (pvs+pas*pva)/det
    omega_am_slow <- (pas*pma)/det
    omega_af_slow <- mu_f*omega_am_slow + (1-mu_f)*omega_au_slow
    omega_ar_slow <- mu_r*omega_am_slow + (1-mu_r)*omega_av_slow
    
    ### fraction of foliage and root litter being transferred to active pool
    transfer_fa <- mu_f*pma + (1-mu_f)*psa
    transfer_ra <- mu_r*pma + (1-mu_r)*psa
    
    ### out df
    ret <- data.frame(decomp_pass, decomp_slow, qq_pass, qq_slow,
                      omega_af_pass, omega_ar_pass, omega_af_slow, omega_ar_slow, 
                      transfer_fa, transfer_ra)
    
    return(ret)
}
