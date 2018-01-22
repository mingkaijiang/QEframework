### NPP as function of nf and LAI (which is calculated from NPP)
### basic function: CUE dependent
eqNC <- function(nf, NPP, CO2, af) {
    
    ##Returns G: total C production (i.e. NPP)
    return(LUE(nf, CO2) * I0 * (1 - exp(-kext*SLA*af*NPP/sf/cfrac)) * cue)
    
}


### NPP as function of nf and LAI (which is calculated from NPP),
### Autotrophic respiration as a function of plant tissue N content
eqNC_respiration <- function(df, NPP, co2, LUE0, Nref, I0, kext, SLA, ADF, sf, cfrac) {
    
    Ra <- Compute_Ra(a=ADF, NPP=NPP)
    
    ##Returns G: total C production (i.e. NPP)
    return(LUE(df, co2, LUE0, Nref) * I0 * (1 - exp(-kext*SLA*ADF$af*NPP/sf/cfrac)) - Ra)
    
}

### NPP as function of nf and LAI (which is calculated from NPP)
### basic function: CUE dependent and based on the full model
### for cn model only
eqPC_full_cn <- function(nf, nfdf, NPP, CO2) {
    
    # in umol C m-2 d-1
    lue_yr <- LUE_full_cn_walker(nf, nfdf, CO2, NPP*1000.0) * par 
    
    # return gpp as kg m-2 yr-1
    gpp <- lue_yr * (1 - exp(-kext*SLA*nfdf$af*NPP/sf/cfrac)) * conv * 365 / 1000.0
    
    ##Returns G: total C production (i.e. NPP)
    return( gpp * cue)
    
}


### NPP as function of nf and LAI (which is calculated from NPP)
### basic function: CUE dependent and based on the full model
### for cn model only
eqPC_simple_cn <- function(nf, nfdf, NPP, CO2) {
    
    # in umol C m-2 d-1
    # lue_yr <- exp((-2.76 + 0.09 * log(nf))) * par
    lue_yr <- exp(-4.23 + 0.25 * log(CO2) + 0.08 * log(nf)) * par
        
    # return gpp as kg m-2 yr-1
    gpp <- lue_yr * (1 - exp(-kext*SLA*nfdf$af*NPP/sf/cfrac)) * conv * 365 / 1000.0
    
    ##Returns G: total C production (i.e. NPP)
    return( gpp * cue)
    
}