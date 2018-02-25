
### This function implements photosynthetic constraint - solve by finding the root
solvePC <- function(nf, af, co2=350,
                    LUE0=1.4, I0=3, Nref=0.04, 
                    kext=0.5, SLA=5, sf=0.5, w = 0.45) {
    # parameters
    # nf is variable
    # making it pass af (fractional allocation to foliage) because this may also be variable
    # co2 = co2 concentration 
    # LUE0 = maximum LUE in kg C GJ-1
    # I0 = total incident radiation in GJ m-2 yr-1
    # Nref = leaf N:C for saturation of photosynthesis
    # kext = light extinction coeffciency
    # SLA = specific leaf area in m2 kg-1 DM
    # sf = turnover rate of foliage in yr-1
    # w = C content of biomass - needed to convert SLA from DM to C
    
    # solve implicit equation
    ans <- c()
    len <- length(nf)
    for (i in 1:len) {
        fPC <- function(NPP) eqPC(nf[i], NPP, co2, LUE0, Nref, I0, kext, SLA, af[i], sf, w) - NPP
        ans[i] <- uniroot(fPC,interval=c(0.1,20))$root
    }
    return(ans)
}