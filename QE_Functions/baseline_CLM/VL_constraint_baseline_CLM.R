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
    NPP_act <- Nmin / soil_n
 
    # return in kg C m-2 yr-1
    NPP <- NPP_act*10^-3     
    
    # out df
    df <- data.frame(NPP,nleach)
    colnames(df) <- c("NPP", "nleach")
    
    return(df)   
}