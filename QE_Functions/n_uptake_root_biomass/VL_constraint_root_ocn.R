
### Calculate the very long term nutrient cycling constraint for N, i.e. passive pool equilibrated
# it is just Nin = Nleach
# specifically for nuptake as a function of root biomass - O-CN approach
# i.e. N uptake as a saturating function of mineral N
VL_constraint_root_ocn <- function(carb_diox) {                   
    # passed are bf and nf, the allocation and plant N:C ratios
    # parameters : 
    # Nin is fixed N inputs (N fixed and deposition) in g m-2 yr-1 (could vary fixation)
    # leachn is the rate of leaching of the mineral N pool (per year)
    # nuptakerate is the rate of N uptake [yr-1] 
    # Nmin is the mineral N pool
    # sr is the decay rate of root in yr-1
    # k - empirically dervied for now
    # vmax
    
    # allocation coefficients
    ar <- aroot
    af <- aleaf
    aw <- 1 - ar - af
    
    # compute Nmin
    Nmin <- Nin / leachn
    
    # arguments
    arg1 <- (ar / sr) * vmax
    arg2 <- Nmin / (Nmin + k)
    arg3 <- arg1 * arg2
    
    # N concentrations of rest of plant   # in g N g-1 C
    if (nwvar == FALSE) {
        nw <- nwood
        nf <- (arg3 - nw * aw) / (af + nrho * ar)
    } else {
        nf <- arg3 / (af + nrho * ar + nwood * aw)
    }
    
    # allocation coefficients
    a_nf <- alloc(nf)
    
    # solve for equilibrium npp
    npp <- photo_constraint_full(nf,a_nf,carb_diox)
    
    ans <- data.frame(nf,npp)
    
    colnames(ans) <- c("nf","NPP")
    
    return(ans)   
}



