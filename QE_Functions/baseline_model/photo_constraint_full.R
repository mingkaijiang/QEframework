### This function implements photosynthetic constraint - solve by finding the root
### Based on the full photosynthesis model for c and n only
photo_constraint_full <- function(nf, nfdf, CO2) {
    # parameters
    # nf is variable
    # making it pass af (fractional allocation to foliage) because this may also be variable
    # CO2 = co2 concentration 
    # kext = light extinction coeffciency
    # SLA = specific leaf area in m2 kg-1 DM
    # sf = turnover rate of foliage in yr-1
    # w = C content of biomass - needed to convert SLA from DM to C
    # cue = carbon use efficiency
    
    ### check length of leaf nc
    len <- length(nf)
    
    ### create output df
    ans <- c()
    
    ### loop to find the unique solution
    for (i in 1:len) {
        fPC <- function(NPP) {
            equil_photo_constraint_full(nf[i], nfdf[i,], NPP, CO2) - NPP
        }
        ans[i] <- uniroot(fPC,interval=c(0.1,20), trace=T)$root
    }
    
    return(ans)
}