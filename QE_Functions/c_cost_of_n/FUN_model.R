FUN_model <- function(nfdf, potnpp) {

    ### plant C:N ratio
    rcn <- 1/(nfdf$nf * nfdf$af + nfdf$nw * nfdf$aw + nfdf$nr * nfdf$ar)

    ### mineral N pool return in kg N m-2 
    Nmin <- Nin / (leachn / (1 - leachn)) *10^-3   

    ### calculate passive N uptake kg N m-2 yr-1
    Npass <- Nmin * (et/sd)
    Ndem <- potnpp / rcn
    Npass2 <- pmin(pmin(Npass, Ndem), Nmin)
    
    ### cost of N fixation, decided not to include this
    # cost_fix <- c_cost_fix()
    
    ### cost of active N uptake kg C m-2
    croot <- potnpp * nfdf$ar * sr 
    cost_active <- c_cost_active(cr=croot)
    
    ### Cost of resorption
    cost_resorb <- c_cost_resorb(potnpp, nfdf)
    
    ### Total C cost
    cost <- pmin(cost_active, cost_resorb)
    
    ### Calculate C spent on uptake N
    Cacq <- (rcn - potnpp / Npass2) / (1 / cost - 1 / Npass2)
    
    ### Calculate real C spent on uptake
    ### incorporating variable CN cost
    modifier <- dynamic_CN_modifier(rcn, cost)
    Cacq_real <- round(Cacq * modifier, 4)
    
    ## NPP allocated to growth
    Cgrow <- round(potnpp,4) - Cacq_real

    ### out df
    out <- data.frame(nfdf$nf, Cgrow, Cacq_real, potnpp, cost)
    colnames(out) <- c("nf", "NPP_grow", "NPP_acq", "NPP_pot", "cost")
    
    return(out)   
}

