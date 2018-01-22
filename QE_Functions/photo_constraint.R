### This function implements photosynthetic constraint - solve by finding the root
photo_constraint <- function(nf, pf, nfdf, pfdf, CO2) {
    # parameters
    # nf is variable
    # making it pass af (fractional allocation to foliage) because this may also be variable
    # co2 = co2 concentration 
    # LUE0 = maximum gross LUE in kg C GJ-1
    # I0 = total incident radiation in GJ m-2 yr-1
    # Nref = leaf N:C for saturation of photosynthesis
    # kext = light extinction coeffciency
    # SLA = specific leaf area in m2 kg-1 DM
    # sf = turnover rate of foliage in yr-1
    # w = C content of biomass - needed to convert SLA from DM to C
    # cue = carbon use efficiency
    
    len <- length(nf)
    
    ans <- c()
    
    for (i in 1:len) {
        fPC <- function(NPP) eqPC(nf[i], pf[i], pfdf$af[i], NPP, CO2) - NPP
        ans[i] <- uniroot(fPC,interval=c(0.1,20), trace=T)$root
        
    }
    
    return(ans)
}

### This function implements photosynthetic constraint - solve by finding the root
### Based on the full photosynthesis model for c and n only
photo_constraint_full_cn <- function(nf, nfdf, CO2) {
    # parameters
    # nf is variable
    # making it pass af (fractional allocation to foliage) because this may also be variable
    # co2 = co2 concentration 
    # LUE0 = maximum gross LUE in kg C GJ-1
    # I0 = total incident radiation in GJ m-2 yr-1
    # Nref = leaf N:C for saturation of photosynthesis
    # kext = light extinction coeffciency
    # SLA = specific leaf area in m2 kg-1 DM
    # sf = turnover rate of foliage in yr-1
    # w = C content of biomass - needed to convert SLA from DM to C
    # cue = carbon use efficiency
    
    len <- length(nf)
    
    ans <- c()
    
    for (i in 1:len) {
        fPC <- function(NPP) eqPC_full_cn(nf[i], nfdf[i,], NPP, CO2) - NPP
        ans[i] <- uniroot(fPC,interval=c(0.1,20), trace=T)$root
    }

    return(ans)
}

### This function implements photosynthetic constraint - solve by finding the root
### Based on the full photosynthesis model for c and n only
photo_constraint_simple_cn <- function(nf, nfdf, CO2) {
    # parameters
    # nf is variable
    # making it pass af (fractional allocation to foliage) because this may also be variable
    # co2 = co2 concentration 
    # LUE0 = maximum gross LUE in kg C GJ-1
    # I0 = total incident radiation in GJ m-2 yr-1
    # Nref = leaf N:C for saturation of photosynthesis
    # kext = light extinction coeffciency
    # SLA = specific leaf area in m2 kg-1 DM
    # sf = turnover rate of foliage in yr-1
    # w = C content of biomass - needed to convert SLA from DM to C
    # cue = carbon use efficiency
    
    len <- length(nf)
    
    ans <- c()
    
    for (i in 1:len) {
        fPC <- function(NPP) eqPC_simple_cn(nf[i], nfdf[i,], NPP, CO2) - NPP
        ans[i] <- uniroot(fPC,interval=c(0.1,20), trace=T)$root
    }
    
    return(ans)
}