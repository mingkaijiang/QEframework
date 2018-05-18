VL_constraint_c_cost <- function(nfdf, potnpp) {

    rcn <- 1/(nfdf$nf * nfdf$af + nfdf$nw * nfdf$aw + nfdf$nr + nfdf$ar)
        
    # mineral N pool return in kg N m-2 
    #Nmin <- Nin  / leachn
    Nmin <- 1
    
    # calculate passive N uptake kg N m-2 yr-1
    Npass <- Nmin * (et/sd)
    Ndem <- potnpp * (1/rcn)
    Npass2 <- pmin(pmin(Npass, Ndem), Nmin)
    
    # cost of N fixation
    cost_fix <- c_cost_fix()
    
    # cost of active N uptake kg C m-2
    croot <- potnpp * nfdf$ar / sr 
    cost_active <- c_cost_active(cr=croot)
    
    # cost of retranslocation
    cost_resorb <- c_cost_resorb(potnpp, nfdf)
    
    # check for the minimum of the three cost
    cost_aq <- pmin(cost_fix, pmin(cost_active, cost_resorb))

    # NPP used for growth
    NPP_grow <- pmax(0, (potnpp - cost_aq * rcn) / (1 - (cost_aq / Npass2)))
    
    # out df
    out <- data.frame(NPP_grow, potnpp, cost_aq)
    colnames(out) <- c("NPPgrow", "NPPpot", "cost")
    
    return(out)   
}