### Allocation and plant N concentrations - required for both PS constraint and NC constraint
allocn <- function(nf,  nwvar = TRUE,
                   nwood = 0.005, #nrho = 0.7,
                   nretrans = 0.5) {
    # parameters
    # nf is the NC ratio of foliage
    # nw is the NC ratio of wood if fixed; otherwise the ratio of wood N:C to foliage N:C
    # nwvar is whether or not to allow wood NC to vary
    # nrho is the ratio of root N:C to foliage N:C
    # nretrans is the fraction of foliage N:C retranslocated  
    
    len <- length(nf)
    ar <- af <- aw <- nw <- nr <- rep(0,len)  # initialise
    for (i in 1:len) {
        ar[i] <- 0.15
        af[i] <- 0.15
        aw[i] <- 1 - ar[i] - af[i]
    }
    
    # N concentrations of rest of plant   # in g N g-1 C
    if (nwvar == FALSE) {
        nw <- nwood
    } else {
        nw <- nwood*nf 
    }
    #nrho <- 0.7
    nr <- nrho*nf
    nfl <- (1.0-nretrans)*nf     
    
    ret <- data.frame(nf,nfl,nw,nr,af,aw,ar)
    return(ret)
}

### Allocation and plant P concentrations - required for both PS constraint and PC constraint
allocp <- function(pf, pwvar) {
    # parameters
    # pf is the PC ratio of foliage
    # pw is the PC ratio of wood if fixed; otherwise the ratio of wood P:C to foliage P:C
    # pwvar is whether or not to allow wood PC to vary
    # prho is the ratio of root P:C to foliage P:C
    # pretrans is the fraction of foliage P:C retranslocated  
    
    len <- length(pf)
    ar <- af <- aw <- pw <- pr <- rep(0,len)  # initialise
    for (i in 1:len) {
        ar[i] <- 0.15
        af[i] <- 0.15
        aw[i] <- 1 - ar[i] - af[i]
    }
    
    # P concentrations of rest of plant   # in g P g-1 C
    if (pwvar == FALSE) {
        pw <- pwood
    } else {
        pw <- pwood*pf 
    }
    #prho <- 0.7
    pr <- prho*pf
    pfl <- (1.0-pretrans)*pf
    
    ret <- data.frame(pf,pfl,pw,pr,af,aw,ar)
    return(ret)
}

### Make inference of pf based on nf
inferpfVL <- function(nf, a) {
    
    # output nf, based on F(nf) = F(pf)
    pf <- c()
    
    Nleach <- (leachn/(1-leachn)) * (a$nfl * a$af + a$nr * a$ar +
                                         a$nw *a$aw)
    
    Pleach <- (leachp/(1-leachp-k1)) 
    Pocc <- (k3/(k2+k3))*(k1/(1-k1-leachp)) 
    
    Pg <- ((Pin * Nleach)/Nin) / (Pocc + Pleach)
    
    if(pwvar == FALSE) {
        pf <- (Pg - pwood * a$aw) / ((1.0 - pretrans) * a$af + prho * a$ar)
    } else {
        pf <- Pg / ((1.0 - pretrans) * a$af + prho * a$ar + a$aw * pwood)
    }
    return(round(pf,8))
}

assim <- function(ci, gamma_star, a1, a2) {
    #    Morning and afternoon calcultion of photosynthesis with the
    #    limitation defined by the variables passed as a1 and a2, i.e. if we
    #     are calculating vcmax or jmax limited.
    
    #    Parameters:
    #    ----------
    #    ci : float
    #    intercellular CO2 concentration.
    #    gamma_star : float
    #    CO2 compensation point in the abscence of mitochondrial respiration
    #    a1 : float
    #    variable depends on whether the calculation is light or rubisco
    #    limited.
    #    a2 : float
    #    variable depends on whether the calculation is light or rubisco
    #    limited.
    
    #    Returns:
    #    -------
    #    assimilation_rate : float
    #    assimilation rate assuming either light or rubisco limitation.
    
    if (ci < gamma_star)
        return (0.0)
    else
        return (a1 * (ci - gamma_star) / (a2 + ci));
    
}


epsilon_simplified <- function(asat, par, alpha, daylen) {
    # simplified function for computing LUE, for using annual results
    lue <- 2.595e-2 + 8.955e-4 * asat - 1.53e-3 * par - 1.118e-1 * alpha
    return(lue)
    
}



### Following two functions calculate NPP - will later need to be replaced by full model
### LUE function of N, P & Ca, based on full photosynthesis model
LUE_full_cnp_walker <- function(nf, pfdf, pf, CO2, NPP) {
    
    ncontent <- NPP * pfdf$af / sf * nf  # g N m-2
    pcontent <- NPP * pfdf$af / sf * pf
    
    # update sla unit from m2 kg-1 DM to m2 g-1
    sla <- SLA / 1000.0
    
    N0 <- ncontent * kn / (1.0 - exp(-kn * sla*pfdf$af*NPP/sf/cfrac))
    P0 <- pcontent * kn / (1.0 - exp(-kn * sla*pfdf$af*NPP/sf/cfrac))
    
    browser()
    #gamma_star <- arrh(mt, gamstar25, eag, tk)
    gamma_star <- 32.97
    
    # Michaelis-Menten coefficents for carboxylation by Rubisco 
    #Kc <- arrh(mt, kc25, eac, tk)
    Kc <- 234.72
    
    # Michaelis-Menten coefficents for oxygenation by Rubisco 
    #Ko <- arrh(mt, ko25, eao, tk)
    Ko <- 216876.747
    
    # return effective Michaelis-Menten coefficient for CO2 
    #km <- (Kc * (1.0 + oi / Ko))
    km <- 461.998
    
    # Walker relationship
    log_vcmax <- 3.946 + 0.921 * log(N0) + 0.121 * log(P0) + 0.282 * log(N0) * log(P0)
    vcmax <- exp(log_vcmax)
    log_jmax <- 1.246 + 0.886 * log_vcmax + 0.089 * log(P0)
    jmax <- exp(log_jmax)
    
    # calculate ci
    g1w <- g1 * wtfac_root
    cica <- g1w / (g1w + sqrt(vpd * PA_2_KPA))
    ci <- cica * CO2
    
    # calculate alpha: quantum efficiency
    alpha <- assim(ci, gamma_star, alpha_j/4.0, 2.0*gamma_star)
    
    ac = assim(ci, gamma_star, vcmax, km)
    
    aj = assim(ci, gamma_star, jmax/4.0, 2.0*gamma_star)
    
    asat <- pmin(aj, ac)
    
    lue_calc <- epsilon_simplified(asat, PAR_MJ, alpha, daylen)
    
    return(lue_calc)
}

cue <- 0.5
kext <- 0.5
nfseq <- seq(0.001, 0.1, by = 0.001)
nrho <- 1.0
CO2 <- 400.0
NPP <- 0.1
PAR_MJ <- 4.0
J_2_UMOL <- 4.57
MJ_TO_J <- 1000000.0
par <- MJ_TO_J * J_2_UMOL * PAR_MJ
sf <- 0.5
SLA <- 5.1
kn <- 0.3
cfrac <- 0.45
oi <- 210000.0
vpd <- 2.4
PA_2_KPA <- 0.001
wtfac_root <- 1.0
g1 <- 3.8667
alpha_j <- 0.308
UMOL_TO_MOL <- 0.000001
MOL_C_TO_GRAMS_C <- 12.0
conv <- UMOL_TO_MOL * MOL_C_TO_GRAMS_C
leachn = 0.05
leachp = 0.1
Nin = 0.2
Pin = 0.02        #0.0086
k1=0.048            # 0.096  # default in full GDAY = 0.048, and 0.4-0.8 of the labile is available to uptake
k2=0.0146           # 0.0146 At EucfACE, based on pH = 4.5, k2 = 0.0146 yr-1, and only half into labile
k3=0.05
a_nf <- as.data.frame(allocn(nfseq,nwvar=FALSE))
pwvar <- FALSE
pwood <- 0.00014
pretrans <- 0.6
prho <- 1.0
pfseq <- inferpfVL(nfseq, a_nf)
a_pf <- as.data.frame(allocp(pfseq,pwvar=FALSE))


### NPP as function of nf and LAI (which is calculated from NPP)
### basic function: CUE dependent and based on the full model
eqPC_full_cnp <- function(nf, pf, pfdf, NPP, CO2) {
    
    # in umol C 
    lue_yr <- LUE_full_cnp_walker(nf, pfdf, pf, CO2, NPP*1000.0/365) * par 
    
    # return gpp as kg m-2 yr-1
    gpp <- lue_yr * (1 - exp(-kext*SLA*pfdf$af*NPP*1000/365/sf/cfrac)) * conv * 365 / 1000.0
    
    ##Returns G: total C production (i.e. NPP)
    return( gpp * cue)
    
}




nfseq <- 0.1
a_nf <- allocn(nfseq, nwvar=FALSE)
pfseq <- inferpfVL(nfseq, a_nf)
a_pf <- as.data.frame(allocp(pfseq,pwvar=FALSE))
NPP <- seq(0.1, 20)

len <- length(nfseq)

ans <- c()

eqPC_full_cnp(nfseq[i], pfseq[i], a_pf[i,], NPP, CO2) - NPP



rm(list=ls(all=TRUE))
