VL_constraint_baseline_CLM <- function(df, a) {
    # passed are bf and nf, the allocation and plant N:C ratios
    # parameters : 
    # Nin is fixed N inputs (N fixed and deposition) in g m-2 yr-1 (could vary fixation)
    # leachn is the rate of leaching of the mineral N pool (per year)
    # nup is the rate of N uptake [yr-1] 
    # Nmin is the mineral N pool
    
    # equation for N constraint with just leaching
    U0 <- Nin
    
    nleach <- leachn
    
    # pool of Nmin
    Nmin <- U0 / nleach
    
    # Calculate NPP
    NPP_NC <- Nmin * nup / (a$nfl*a$af + a$nr*a$ar + a$nw*a$aw)
    # The above equation is the same as:
    # NPP_NC <- U0 * nuptakerate / (nleach * (a$nfl*a$af + a$nr*(a$ar) + a$nw*a$aw))
    
    # return in kg C m-2 yr-1
    NPP_N <- NPP_NC*10^-3     
    
    # out df
    df <- data.frame(NPP_N,nleach)
    colnames(df) <- c("NPP", "nleach")
    
    return(df)   
}