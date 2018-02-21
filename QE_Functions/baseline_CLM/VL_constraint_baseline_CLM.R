VL_constraint_baseline_CLM <- function(df,a) {
    # passed are bf and nf, the allocation and plant N:C ratios
    # parameters : 
    # Nin is fixed N inputs (N fixed and deposition) in g m-2 yr-1 (could vary fixation)
    # leachn is the rate of leaching of the mineral N pool (per year)
    # nup is the rate of N uptake [yr-1] 
    # Nmin is the mineral N pool
    
    # equation for N constraint with just leaching
    U0 <- Nin
    
    # Calculate NPP
    nleach <- leachn/(1-leachn) * (a$nfl*a$af + a$nr*a$ar + a$nw*a$aw)
 
    # return in kg C m-2 yr-1
    NPP <- U0 / nleach *10^-3     
    
    # out df
    out <- data.frame(NPP,nleach)
    colnames(out) <- c("NPP", "nleach")
    
    return(out)   
}