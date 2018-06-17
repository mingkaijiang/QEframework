FUN_model_1 <- function(nfdf, potnpp) {

    ### plant C:N ratio
    #rcn <- 1/(nfdf$nf * nfdf$af + nfdf$nw * nfdf$aw + nfdf$nr * nfdf$ar)
    rcn <- 1/nfdf$nf
    
    ### mineral N pool return in kg N m-2 
    Nmin <- Nin *10^-3 / (leachn / (1 - leachn)) 

    ### calculate passive N uptake kg N m-2 yr-1
    Npass <- Nmin * (et/sd)
    Ndem <- potnpp / rcn
    Npass2 <- pmin(pmin(Npass, Ndem), Nmin)
    
    ### cost of N fixation, decided not to include this
    # cost_fix <- c_cost_fix()
    
    ### cost of active N uptake kg C m-2
    # croot <- potnpp * nfdf$ar / sr 
    # cost_active <- c_cost_active(cr=croot)
    
    ### Cost of resorption
    # cost_resorb <- c_cost_resorb(potnpp, nfdf)
    
    ### Total C cost
    # cost <- pmin(cost_active, cost_resorb)
    cost <- 0.0
    
    ### Calculate C spent on uptake N
    # Cacq <- (potnpp - rcn * Npass2) / (rcn / cost + 1.0)
    Cacq <- 0
    
    ### Calculate real C spent on uptake
    ### incorporating variable CN cost
    # modifier <- dynamic_CN_modifier(rcn, cost)
    # Cacq_real <- round(Cacq * modifier, 4)
    
    ## NPP allocated to growth
    Cgrow <- rcn * Npass2

    ### out df
    out <- data.frame(nfdf$nf, Cgrow, Cacq, potnpp, cost)
    colnames(out) <- c("nf", "NPP_grow", "NPP_acq", "NPP_pot", "cost")
    
    return(out)   
}

