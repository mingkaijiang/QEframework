
### NPP as function of nf and LAI (which is calculated from NPP)
eqPC <- function(nf, NPP, co2, LUE0, Nref, I0, kext, SLA, af, sf, cfrac) {
    
    return(LUE(nf, co2, LUE0, Nref) * I0 * (1 - exp(-kext*SLA*af*NPP/sf/cfrac)))
    
}