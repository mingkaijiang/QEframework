VL_constraint_c_cost <- function(nfdf, potnpp) {

    # mineral N pool return in kg N m-2 
    Nmin <- Nin / 1000 / leachn
    
    # calculate passive N uptake kg N m-2 yr-1
    Npass <- Nmin * (et/sd)
    Ndem <- potnpp * (nfdf$nf * nfdf$af + nfdf$nw * nfdf$aw + nfdf$nr + nfdf$ar)
    Npass2 <- min(Npass, Ndem)
    
    # calculate c acqured due to passive N uptake
    c_pass <- Npass2 / (nfdf$nf * nfdf$af + nfdf$nw * nfdf$aw + nfdf$nr + nfdf$ar)
    
    # cost of N fixation
    cost_fix <- c_cost_fix()
    
    # cost of active N uptake kg C m-2
    croot <- potnpp * nfdf$ar / sr 
    cost_active <- c_cost_active(cr=croot)
    
    # cost of retranslocation
    cost_resorb <- c_cost_resorb(nfdf)
    
    # check for the minimum of the three cost
    cost_aq <- min(cost_fix, min(cost_active, cost_resorb))
    
    # NPP used for growth
    NPP_act <-  potnpp - cost_aq + c_pass
    
    browser()
    
    # out df
    out <- data.frame(NPP_act, cost_aq)
    colnames(out) <- c("NPPact", "NPPaqu")
    
    return(out)   
}