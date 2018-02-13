
### Calculate the very long term nutrient cycling constraint for N, i.e. passive pool equilibrated
# it is just Nin = Nleach
# specifically for nuptake as a function of root biomass - CLM approach
# i.e. N uptake as a saturating function of mineral N
VL_constraint_root_clm <- function(df, a) {                   
    # passed are bf and nf, the allocation and plant N:C ratios
    # parameters : 
    # Nin is fixed N inputs (N fixed and deposition) in g m-2 yr-1 (could vary fixation)
    # leachn is the rate of leaching of the mineral N pool (per year)
    # nuptakerate is the rate of N uptake [yr-1] 
    # Nmin is the mineral N pool
    
    # compute Nmin
    Nmin <- Nin / leachn
    
    # leaf CN scalar
    scalar_n <- (1/df - cn.leaf.min) / (cn.leaf.max - cn.leaf.min)
    
    # arguments
    arg1 <- (a$ar / sr) * umax
    arg2 <- Nmin / (Nmin + ksmin)
    
    # total n uptake as g N m-2 yr-1
    n_uptake <- arg1 * arg2 * scalar_temp * scalar_n
    
    # NPP as g C m-2 yr-1
    NPP <- n_uptake / (a$nfl*a$af + a$nr*a$ar + a$nw*a$aw)
    
    ans <- data.frame(df, NPP)
    
    colnames(ans) <- c("nf", "NPP")
    
    return(ans)   
}



