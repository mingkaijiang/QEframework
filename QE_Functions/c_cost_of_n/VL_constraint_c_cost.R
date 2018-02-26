VL_constraint_c_cost <- function(nfdf, potnpp) {

    # mineral N pool return in kg N m-2 yr-1
    Nmin <- Nin / 1000 / leachn
    
    # cost of N fixation
    cost_fix <- c_cost_fix()
    
    # croot

    # Calculate NPP used for growth
    NPP_growth <- U0 / nleach
    
    # cost of active N uptake
    cost_active <- 0 #(kn2 / Nmin) * (kc2 / croot)
    
    # cost of resorb
    cost_resorb <- kr / nfdf$nf
    
    # Cost of acquiring (kg C kg N-1)
    # minimum of cost of resorption, active uptake or fixation
    cost_aq <- min(cost_resorb, cost_active, cost_fix)
    
    # NPP used for N uptake
    # N acquired * cost of aquiring (kg C kg N-1)
    NPP_aq <- N_aq * cost_aq

    # total potential NPP is the sum of NPP used for growth and NPP used for N acquisition
    NPP_total <- NPP_growth + NPP_aq
    
    # return in kg C m-2 yr1
    NPP <- NPP_total * 10^-3
    
    # out df
    df <- data.frame(NPP, NPP_growth, NPP_aq, nleach)
    colnames(df) <- c("NPP", "NPP_growth", "NPP_aq", "nleach")
    
    return(df)   
}