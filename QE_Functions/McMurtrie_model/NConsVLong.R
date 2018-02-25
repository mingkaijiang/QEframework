### Calculate the very long term nutrient cycling constraint, i.e. passive pool equilibrated
# it is just Nin = Nleach
NConsVLong <- function(nf, a, Nin=1.0, leach=0.05) {
    # passed are nf and a, the allocation and plant N:C ratios
    # parameters : 
    # Nin is fixed N inputs (N deposition, N fixation) in g m-2 yr-1 (could vary fixation)
    # leach is the rate of leaching of the mineral pool (per year)
    
    # equation for N constraint with just leaching
    U0 <- Nin
    nleach <- leach/(1-leach) * (a$nfl*a$af + a$nr*(a$ar) + a$nw*a$aw)
    NPP_NC <- U0 / (nleach)   # will be in g C m-2 yr-1
    NPP <- NPP_NC*10^-3 # returned in kg C m-2 yr-1
    df <- data.frame(NPP,nleach)
    return(df)   
}
