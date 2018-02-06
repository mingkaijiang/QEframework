
# Find the medium term equilibrium nf and NPP under standard conditions - by finding the root
# specifically for original model (i.e. Figure 1)
solveMedium_original <- function(CO2,Cpass,Cslow,NinL) {
    fn <- function(nf) {
        solveNC(nf,allocn(nf)$af,CO2) - NConsMedium(nf,allocn(nf),Cpass,Cslow,NinL)$NPP
        
    }
    equilnf <- uniroot(fn,interval=c(0.01,0.05))$root
    equilNPP <- solveNC(equilnf,af=allocn(equilnf)$af, CO2)
    
    ans <- data.frame(equilnf,equilNPP)
    return(ans)
}


# Find the medium term equilibrium nf and NPP under standard conditions - by finding the root
# specifically for exudation model, without considering exudation
solveMedium <- function(CO2,Cpass,Cslow,NinL) {
    fn <- function(nf) {
        photo_constraint_full_cn(nf, allocn(nf),
                                 CO2) - NConsMedium(nf,allocn(nf),Cpass,Cslow,NinL)$NPP
        
    }
    equilnf <- uniroot(fn,interval=c(0.001,0.1))$root
    equilNPP <- photo_constraint_full_cn(equilnf, 
                                         allocn(equilnf), CO2)
    # browser()
    ans <- data.frame(equilnf,equilNPP)
    return(ans)
}


# Find the medium term equilibrium nf and NPP under standard conditions - by finding the root
# specifically for exudation model, turning exudation on
solveMedium_priming <- function(CO2,Cpass,Cslow,NinL) {
    fn <- function(nf) {
        photo_constraint_full_cn(nf, allocn_exudation(nf),
                                 CO2) - NConsMedium_priming(nf,allocn_exudation(nf),Cpass,Cslow,NinL)$NPP
        
    }
    equilnf <- uniroot(fn,interval=c(0.001,0.1))$root
    equilNPP <- photo_constraint_full_cn(equilnf, 
                                         allocn_exudation(equilnf), CO2)
    ans <- data.frame(equilnf,equilNPP)
    return(ans)
}


# Find the medium term equilibrium nf and NPP under standard conditions - by finding the root
# specifically for explicit mineral pools
solveMedium_expl_min <- function(CO2, Cpass, Cslow, NinL) {
    fn <- function(nf) {
        photo_constraint_full_cn(nf, 
                                  allocn(nf), 
                                  CO2) - NConsMedium_expl_min(nf,allocn(nf),Cpass,Cslow, NinL)$NPP
        
    }
    equilnf <- uniroot(fn,interval=c(0.001,0.1))$root
    equilNPP <- photo_constraint_full_cn(equilnf, 
                                          allocn(equilnf),  CO2)
    ans <- data.frame(equilnf,equilNPP)
    return(ans)
}

# Find the medium term equilibrium nf and NPP under standard conditions - by finding the root
# specifically for root gday model
solveMedium_root_gday <- function(CO2,Cpass,Cslow,NinL) {
    fn <- function(nf) {
        photo_constraint_full_cn(nf, allocn(nf),
                                 CO2) - NConsMedium_root_gday(nf,allocn(nf),NinL,Cpass,Cslow)$NPP
        
    }
    equilnf <- uniroot(fn,interval=c(0.001,0.1))$root
    equilNPP <- photo_constraint_full_cn(equilnf, 
                                         allocn(equilnf), CO2)
    ans <- data.frame(equilnf,equilNPP)
    return(ans)
}

# Find the medium term equilibrium nf and NPP under standard conditions - by finding the root
# specifically for root ocn model
solveMedium_root_ocn <- function(CO2,Cpass,Cslow,NinL) {
    fn <- function(nf) {
        photo_constraint_full_cn(nf, allocn(nf),
                                 CO2) - NConsMedium_root_ocn(nf,allocn(nf),Cpass,Cslow,NinL)$NPP
        
    }
    equilnf <- uniroot(fn,interval=c(0.001,0.1))$root
    equilNPP <- photo_constraint_full_cn(equilnf, 
                                         allocn(equilnf), CO2)
    ans <- data.frame(equilnf,equilNPP)
    return(ans)
}