### Compute CUE at various timesteps
cue_compute <- function(nf, pf, nfdf, pfdf, NPP, CO2) {
    
    lue_yr <-  LUE_full_cnp_walker(nf, pfdf, pf, CO2, NPP*1000.0) * par
    
    Rd <- Compute_Rdark(nfdf, pfdf, NPP*1000.0)
    
    Ra <- Rd * 2.5
    
    # return gpp as kg m-2 yr-1
    gpp <- lue_yr * (1 - exp(-kext*SLA*nfdf$af*NPP/sf/cfrac)) * conv * 365 / 1000.0
    
    ##Returns G: total C production (i.e. NPP)
    #return(NPP/gpp)
    return((gpp-Ra)/gpp)
}