c_cost_active <- function(cr) {
    ### Calculate c cost of active N uptake
    ### parameters:
    ### kn2 and kc2 are both 1 kg C m-2, derivations available in Fisher et al. 2010
    ### c_cost_active: kg C kg N-1
    
    Nmin <- Nin / 1000 / leachn
    
    c_cost_active <- (kn2 / Nmin) * (kc2 / cr)
    
    return(c_cost_active)
}