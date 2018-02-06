# Find the long term equilibrium nf and NPP under standard conditions - by finding the root
solveLong <- function(CO2,Cpass,NinL, PinL) {
    fn <- function(nf) {
        photo_constraint(nf, inferpfVL(nf, allocn(nf)), 
                         allocn(nf), 
                         allocp(inferpfVL(nf, allocn(nf))), 
                         CO2) - Long_constraint_N(nf,allocn(nf),Cpass=Cpass,NinL)$NPP
    }
    equilnf <- uniroot(fn,interval=c(0.01,0.1))$root
    equilpf <- inferpfVL(equilnf, allocn(equilnf))
    equilNPP <- photo_constraint(equilnf, equilpf, 
                                 allocn(equilnf), allocp(equilpf), CO2)
    ans <- data.frame(equilnf, equilpf, equilNPP)
    return(ans)
}

# Find the long term equilibrium nf and NPP under standard conditions - by finding the root
# specifically for considering both long and medium terms together
solveLong_full_cn_medium <- function(CO2,Cpass,NinL) {
    fn <- function(nf) {
        photo_constraint_full_cn(nf, allocn_exudation(nf), 
                                 CO2) - Long_constraint_N(nf,allocn_exudation(nf),Cpass=Cpass,NinL)$NPP
    }
    equilnf <- uniroot(fn,interval=c(0.001,0.1))$root
    equilNPP <- photo_constraint_full_cn(equilnf, 
                                         allocn_exudation(equilnf), CO2)
    equilpf <- "NA"
    ans <- data.frame(equilnf, equilpf, equilNPP)
    return(ans)
}

# Find the long term equilibrium nf and NPP under standard conditions - by finding the root
solveLong_simple_cn <- function(CO2,Cpass,NinL) {
    fn <- function(nf) {
        photo_constraint_simple_cn(nf, allocn(nf), 
                                 CO2) - Long_constraint_N(nf,allocn(nf),Cpass=Cpass,NinL)$NPP
    }
    equilnf <- uniroot(fn,interval=c(0.01,0.1))$root
    equilNPP <- photo_constraint_simple_cn(equilnf, 
                                         allocn(equilnf), CO2)
    equilpf <- "NA"
    ans <- data.frame(equilnf, equilpf, equilNPP)
    return(ans)
}


solveLong_expl_min <- function(CO2,Cpass,NinL) {
    fn <- function(nf) {
        photo_constraint_full_cn(nf, 
                                  allocn(nf), 
                                  CO2) - NConsLong_expl_min(nf,allocn(nf),Cpass,NinL)$NPP
        
    }
    equilnf <- uniroot(fn,interval=c(0.001,0.1))$root
    equilNPP <- photo_constraint_full_cn(equilnf, 
                                          allocn(equilnf),  CO2)
    ans <- data.frame(equilnf,equilNPP)
    return(ans)
}

# Find the long term equilibrium nf and NPP under standard conditions - by finding the root
solveLong_respiration <- function(CO2,Cpass,NinL, PinL) {
    fn <- function(nf) {
        photo_constraint_respiration(nf, inferpfVL(nf, allocn(nf)), 
                                     allocn(nf), 
                                     allocp(inferpfVL(nf, allocn(nf))), 
                                     CO2) - Long_constraint_N(nf,allocn(nf),Cpass=Cpass,NinL)$NPP
    }
    equilnf <- uniroot(fn,interval=c(0.0085,0.04))$root
    equilpf <- inferpfVL(equilnf, allocn(equilnf))
    equilNPP <- photo_constraint_respiration(equilnf, equilpf, 
                                          allocn(equilnf), allocp(equilpf), CO2)
    ans <- data.frame(equilnf, equilpf, equilNPP)
    return(ans)
}

# Find the long term equilibrium nf and NPP under standard conditions - by finding the root
# specifically for nuptake ~ root biomass  - O-CN approach
# i.e. N uptake as a saturating function of mineral N
solveLong_root_ocn <- function(CO2, Cpass, NinL) {
    fn <- function(nf) {
        photo_constraint_full_cn(nf, allocn(nf), CO2) - NConsLong_root_ocn(nf,allocn(nf),Cpass=Cpass,NinL)$NPP
    }
    equilnf <- uniroot(fn,interval=c(0.01,0.1))$root
    equilNPP <- photo_constraint_full_cn(equilnf, allocn(equilnf), CO2)
    ans <- data.frame(equilnf,equilNPP)
    return(ans)
}


# Find the long term equilibrium nf and NPP under standard conditions - by finding the root
# specifically for nuptake ~ root biomass  - O-CN approach
# i.e. N uptake as a saturating function of mineral N
solveLong_root_gday <- function(CO2,Cpass,NinL) {
    fn <- function(nf) {
        photo_constraint_full_cn(nf, allocn(nf), CO2) - NConsLong_root_gday(nf,allocn(nf),Cpass=Cpass,NinL)$NPP
    }
    equilnf <- uniroot(fn,interval=c(0.001,0.1))$root
    equilNPP <- photo_constraint_full_cn(equilnf, 
                                          allocn(equilnf), CO2)
    ans <- data.frame(equilnf,equilNPP)
    return(ans)
}

# Find the long term equilibrium nf and NPP under standard conditions - by finding the root
# specifically for exudation model
solveLong_exudation <- function(CO2,Cpass,NinL) {
    fn <- function(nf) {
        photo_constraint_full_cn(nf, allocn_exudation(nf),
                                 CO2) - NConsLong_exudation(nf,allocn_exudation(nf),Cpass,NinL)$NPP
        
    }
    equilnf <- uniroot(fn,interval=c(0.004,0.1))$root
    equilNPP <- photo_constraint_full_cn(equilnf, 
                                          allocn_exudation(equilnf), CO2)
    ans <- data.frame(equilnf,equilNPP)
    return(ans)
}

# Find the long term equilibrium nf and NPP under standard conditions - by finding the root
# specifically for exudation model with medium pool considered
solveLong_exudation_medium <- function(CO2,Cpass,NinL) {
    fn <- function(nf) {
        photo_constraint_full_cn(nf, allocn_exudation(nf),
                                 CO2) - NConsLong_exudation_medium(nf,allocn_exudation(nf),Cpass,NinL)$NPP
        
    }
    equilnf <- uniroot(fn,interval=c(0.004,0.1))$root
    equilNPP <- photo_constraint_full_cn(equilnf, 
                                         allocn_exudation(equilnf), CO2)
    ans <- data.frame(equilnf,equilNPP)
    return(ans)
}