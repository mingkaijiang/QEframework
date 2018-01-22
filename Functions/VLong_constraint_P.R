### Calculate the very long term nutrient cycling constraint for P, i.e. passive pool equilibrated
# it is just Pin = Pleach + Pocc
VLong_constraint_P <- function(pf, pfdf) {
    # parameters : 
    # Pin is P deposition inputs in g m-2 yr-1 (could vary fixation)
    # leachp is the rate of leaching of the labile P pool (per year)
    # k1 is the transfer rate from labile to secondary P pool
    # k2 is the transfer rate from secondary to labile P pool
    # k3 is the transfer rate from secondary to occluded P pool
    
    U0 = Pin
    pleach <- (leachp/(1-leachp-k1)) * (pfdf$pfl*pfdf$af + pfdf$pr*pfdf$ar + pfdf$pw*pfdf$aw)
    pocc <- (k3/(k2+k3))*(k1/(1-k1-leachp)) * (pfdf$pfl*pfdf$af + pfdf$pr*pfdf$ar + pfdf$pw*pfdf$aw)
    
    NPP_PC <- U0 / (pleach + pocc)   # will be in g C m-2 yr-1
    NPP_P <- NPP_PC*10^-3     # returned in kg C m-2 yr-1
    
    df <- data.frame(NPP_P,pleach, pocc)
    return(df)   
}