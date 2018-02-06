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


