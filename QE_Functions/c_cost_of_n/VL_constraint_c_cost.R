VL_constraint_c_cost <- function(nfdf, potnpp) {

    ### plant C:N ratio
    rcn <- 1/(nfdf$nf * nfdf$af + nfdf$nw * nfdf$aw + nfdf$nr * nfdf$ar)
    #rcn <- seq(300, 10, length.out = length(potnpp))
    
    ### mineral N pool return in kg N m-2 
    Nmin <- Nin / leachn

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
    
    ## NPP allocated to growth
    Cgrow <- potnpp - Cacq
  
    ### out df
    out <- data.frame(Cgrow, Cacq, potnpp, cost)
    colnames(out) <- c("NPP_grow", "NPP_acq", "NPP_pot", "cost")
    
    return(out)   
}

### a summary of issues:
### 3. Variable plant N:C ratio not considered. Need to incorporate this. 
