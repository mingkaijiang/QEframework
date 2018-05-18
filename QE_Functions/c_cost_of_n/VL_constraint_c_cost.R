VL_constraint_c_cost <- function(nfdf, potnpp) {

    rcn <- 1/(nfdf$nf * nfdf$af + nfdf$nw * nfdf$aw + nfdf$nr * nfdf$ar)
    # rcn <- 300
    
    # mineral N pool return in kg N m-2 
    Nmin <- Nin
    #Nmin <- 1
    
    # calculate passive N uptake kg N m-2 yr-1
    Npass <- Nmin * (et/sd)
    Ndem <- potnpp / rcn
    Npass2 <- pmin(pmin(Npass, Ndem), Nmin)
    
    # cost of N fixation
    cost_fix <- c_cost_fix()
    
    # cost of active N uptake kg C m-2
    croot <- potnpp * nfdf$ar * sr 
    # croot <- 0.2
    cost_active <- c_cost_active(cr=croot)
    
    # N conductance  kgN/kgC
    Nconduct <- sum(1/cost_fix, 1/cost_active)
    
    # fraction of C allocated to each pathway
    Cfrac_fix <- (1/cost_fix)/Nconduct
    Cfrac_active <- (1/cost_active)/Nconduct
    
    # N acquired from each stream per unit C spent (kgN/kgC):
    Nextra_fix <- Cfrac_fix / cost_fix
    Nextra_active <- Cfrac_active / cost_active
    
    # total amount of N uptake per unit C spent (gN/gC):
    Nextra_tot <- Nextra_fix + Nextra_active
    
    # overal N cost is:
    Ncost_tot <- 1/Nextra_tot
    
    ## NPP_grow
    NPP_grow <- ((potnpp) * rcn) / (rcn + Ncost_tot) 
    
    # NPP acquisition
    NPP_aq <- potnpp - NPP_grow
  
    # out df
    out <- data.frame(NPP_grow, NPP_aq, potnpp, cost_aq)
    colnames(out) <- c("NPPgrow", "NPP_acq", "NPPpot", "cost")
    
    return(out)   
}

### a summary of issues:
### 1. Nmin is not fully considered. Right now I am using a simplified Nmin of Nmin = Nin. So leaching was not considered. 
### 2. Effect of passive uptake not considered. Passive uptake means that plants are taking up N for free, hence there should have
###    more NPP allocated to growth, but right now the implementation doesn't include free N.
### 3. Leaf N retranslocation is not yet considered. The online FUN model descrption has two parts: free and paid-for N retranslocation. 
###    It does seem quite challenging to implement, but is not impossible. Can refer to Fisher et al. 2010 for a simplified version. 
###    But unit there is incorrect. 
