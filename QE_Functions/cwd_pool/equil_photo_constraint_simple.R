### NPP as function of nf and LAI (which is calculated from NPP)
### basic function: CUE dependent and based on the simple model
### for cn model only
equil_photo_constraint_simple <- function(nf, nfdf, NPP, CO2) {
    
    # in umol C m-2 d-1
    # lue_yr <- exp((-2.76 + 0.09 * log(nf))) * par
    lue_yr <- exp(-4.23 + 0.25 * log(CO2) + 0.08 * log(nf)) * par
    
    # return gpp as kg m-2 yr-1
    gpp <- lue_yr * (1 - exp(-kext*SLA*nfdf$af*NPP/sf/cfrac)) * conv * 365 / 1000.0
    
    ##Returns G: total C production (i.e. NPP)
    return( gpp * cue)
    
}