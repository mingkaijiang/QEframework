### This function implements photosynthetic constraint - solve by finding the root
### Based on the FUN photosynthesis model for c and n only
photo_constraint_FUN <- function(nf, nfdf, CO2) {
    
    ### check length of leaf nc
    len <- length(nf)
    
    ### create output df
    ans <- c()
    
    ### loop to find the unique solution
    for (i in 1:len) {
        fPC <- function(NPP) {
            equil_photo_constraint_full(nf[i], nfdf[i,], NPP, CO2) - NPP
        }
        ans[i] <- uniroot(fPC,interval=c(0.1,20), trace=T)$root
    }
    
    potnpp <- ans
    
    ### plant C:N ratio
    rcn <- 1/(nfdf$nf * nfdf$af + nfdf$nw * nfdf$aw + nfdf$nr * nfdf$ar)
    
    ### mineral N pool return in kg N m-2 
    Nmin <- Nin *10^-3 / (leachn / (1 - leachn)) 
    
    ### calculate passive N uptake kg N m-2 yr-1
    Npass <- Nmin * (et/sd)
    Ndem <- potnpp / rcn
    Npass2 <- pmin(pmin(Npass, Ndem), Nmin)
    
    ### cost of N fixation, decided not to include this
    cost_fix <- c_cost_fix()
    
    ### cost of active N uptake kg C m-2
    croot <- potnpp * nfdf$ar / sr 
    cost_active <- c_cost_active(cr=croot)
    
    ### Cost of resorption
    cost_resorb <- c_cost_resorb(potnpp, nfdf)
    
    ### Total C cost
    cost <- pmin(cost_active, pmin(cost_resorb, cost_fix))
    
    ### Calculate C spent on uptake N
    Cacq <- (potnpp - rcn * Npass2) / (rcn / cost + 1.0)
    
    ### Calculate real C spent on uptake
    ### incorporating variable CN cost
    modifier <- dynamic_CN_modifier(rcn, cost)
    Cacq_real <- round(Cacq * modifier, 4)
    
    ## NPP allocated to growth
    Cgrow <- round(potnpp,4) - Cacq_real
    
    ### out df
    out <- data.frame(nfdf$nf, Cgrow, Cacq_real, potnpp, cost, cost_resorb, cost_active, cost_fix)
    colnames(out) <- c("nf", "NPP_grow", "NPP_acq", "NPP_pot", "cost", "cost_resorb", "cost_active", "cost_fix")
    
    return(out)  
    
}