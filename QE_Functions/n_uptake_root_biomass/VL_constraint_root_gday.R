
VL_constraint_root_gday <- function(df, a) {                   
    ### Calculate the very long term nutrient cycling constraint for N, i.e. passive pool equilibrated
    # it is just Nin = Nleach
    # specifically for nuptake as a function of root biomass - GDAY approach
    # i.e. N uptake as a saturating function of root biomass
    # passed are bf and nf, the allocation and plant N:C ratios
    # parameters : 
    # Nin is fixed N inputs (N fixed and deposition) in g m-2 yr-1 (could vary fixation)
    # leachn is the rate of leaching of the mineral N pool (per year)
    # Nmin is the mineral N pool
    # sr is the decay rate of root in yr-1
    # kr is the value of root carbon at which 50% of the available N is taken up  
    
    
    # equation for N constraint with just leaching
    U0 <- Nin
    
    # N min pool
    Nmin <- Nin / leachn
    
    # allocation stuffs
    a_plant <- a$nfl*a$af + a$nr*a$ar + a$nw*a$aw
    
    # root biomass
    root_biomass <- a$ar / sr
    
    # leaching term
    nleach <- Nmin * leachn
    
    # equation for NPP
    #NPP_NC <- (root_biomass * Nmin - (A_NF * kr)) / (A_NF * root_biomass)
    NPP_NC <- (Nmin - kr * a_plant / root_biomass ) / a_plant
    
    # returned in kg C m-2 yr-1
    NPP <- NPP_NC*10^-3     
    
    # out df
    out <- data.frame(NPP, nleach)
    
    return(out)   
}