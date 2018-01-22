# Find the very-long term equilibrium nf and NPP under standard conditions - by finding the root
solveVLong <- function(CO2) {
    fn <- function(nf) {
        photo_constraint(nf, inferpfVL(nf, allocn(nf)), 
                         allocn(nf),allocp(inferpfVL(nf, allocn(nf))), 
                         CO2) - VLong_constraint_N(nf,allocn(nf))$NPP
    }
    equilnf <- uniroot(fn,interval=c(0.01,0.1))$root
    equilpf <- inferpfVL(equilnf, allocn(equilnf))
    equilNPP <- photo_constraint(equilnf, equilpf, 
                                 allocn(equilnf), allocp(equilpf), CO2)
    
    ans <- data.frame(equilnf, equilpf, equilNPP)
    return(ans)
}


# Find the very-long term equilibrium nf and NPP under standard conditions - by finding the root
solveVLong_full_cnp <- function(CO2) {
    fn <- function(nf) {
        photo_constraint_full_cnp(nf, inferpfVL(nf, allocn(nf)), 
                                  allocn(nf),allocp(inferpfVL(nf, allocn(nf))), 
                                  CO2) - VLong_constraint_N(nf,allocn(nf))$NPP
    }
    equilnf <- uniroot(fn,interval=c(0.001,0.1))$root
    equilpf <- inferpfVL(equilnf, allocn(equilnf))
    equilNPP <- photo_constraint_full_cnp(equilnf, equilpf, 
                                          allocn(equilnf), allocp(equilpf), CO2)
    
    ans <- data.frame(equilnf, equilpf, equilNPP)
    return(ans)
}

# Find the very-long term equilibrium nf and NPP under standard conditions - by finding the root
solveVLong_simple_cnp <- function(CO2) {
    fn <- function(nf) {
        photo_constraint_simple_cnp(nf, inferpfVL(nf, allocn(nf)), 
                                  allocn(nf),allocp(inferpfVL(nf, allocn(nf))), 
                                  CO2) - VLong_constraint_N(nf,allocn(nf))$NPP
    }
    equilnf <- uniroot(fn,interval=c(0.01,0.1))$root
    equilpf <- inferpfVL(equilnf, allocn(equilnf))
    equilNPP <- photo_constraint_simple_cnp(equilnf, equilpf, 
                                          allocn(equilnf), allocp(equilpf), CO2)
    
    ans <- data.frame(equilnf, equilpf, equilNPP)
    return(ans)
}

# Find the very-long term equilibrium nf and NPP under standard conditions - by finding the root
# specifically for explicit mineral pools
#solveVLong_expl_min <- function(CO2) {
#    fn <- function(nf) {
#        photo_constraint_full_cnp(nf, inferpfVL_expl_min(nf, allocn(nf)), 
#                                  allocn(nf),allocp(inferpfVL_expl_min(nf, allocn(nf))), 
#                                  CO2) - NConsVLong_expl_min(nf,allocn(nf))$NPP
#    }
#    equilnf <- uniroot(fn,interval=c(0.01,0.1))$root
#    equilpf <- inferpfVL_expl_min(equilnf, allocn(equilnf))
#    equilNPP_N <- photo_constraint_full_cnp(equilnf, equilpf, 
#                                            allocn(equilnf), allocp(equilpf), CO2)
#    
#    ans <- data.frame(equilnf,equilpf,equilNPP_N)
#    return(ans)
#}

solveVLong_expl_min <- function(CO2) {
    fn <- function(nf) {
        photo_constraint_full_cn(nf, 
                                  allocn(nf), 
                                  CO2) - NConsVLong_expl_min(nf,allocn(nf))$NPP
    }
    equilnf <- uniroot(fn,interval=c(0.001,0.1))$root
    equilNPP_N <- photo_constraint_full_cn(equilnf,
                                            allocn(equilnf),  CO2)
    
    ans <- data.frame(equilnf,equilNPP_N)
    return(ans)
}

# Find the very-long term equilibrium nf and NPP under standard conditions - by finding the root
solveVLong_full_cnp_fix_wood <- function(CO2) {
    fn <- function(nf) {
        photo_constraint_full_cnp(nf, inferpfVL(nf, allocn(nf)), 
                                  allocn(nf),allocp(inferpfVL(nf, allocn(nf))), 
                                  CO2) - VLong_constraint_N(nf,allocn(nf))$NPP
    }
    equilnf <- uniroot(fn,interval=c(0.01,0.1))$root
    equilpf <- inferpfVL(equilnf, allocn(equilnf))
    equilNPP <- photo_constraint_full_cnp(equilnf, equilpf, 
                                          allocn(equilnf), allocp(equilpf), CO2)
    
    ans <- data.frame(equilnf, equilpf, equilNPP)
    return(ans)
}

# Find the very-long term equilibrium nf and NPP under standard conditions - by finding the root
solveVLong_full_cn <- function(CO2) {
    fn <- function(nf) {
        photo_constraint_full_cn(nf, allocn(nf),CO2) - VLong_constraint_N(nf,allocn(nf))$NPP
    }
    equilnf <- uniroot(fn,interval=c(0.001,0.1))$root
    equilNPP <- photo_constraint_full_cn(equilnf, allocn(equilnf), CO2)
    equilpf <- "NA"
    ans <- data.frame(equilnf, "NA", equilNPP)
    return(ans)
}

# Find the very-long term equilibrium nf and NPP under standard conditions - by finding the root
solveVLong_simple_cn <- function(CO2) {
    fn <- function(nf) {
        photo_constraint_simple_cn(nf, allocn(nf),CO2) - VLong_constraint_N(nf,allocn(nf))$NPP
    }
    equilnf <- uniroot(fn,interval=c(0.01,0.1))$root
    equilNPP <- photo_constraint_simple_cn(equilnf, allocn(equilnf), CO2)
    equilpf <- "NA"
    ans <- data.frame(equilnf, "NA", equilNPP)
    return(ans)
}

# Find the very-long term equilibrium nf and NPP under standard conditions - by finding the root
solveVLong_respiration <- function(CO2) {
    fn <- function(nf) {
        photo_constraint_respiration(nf, inferpfVL(nf, allocn(nf)), 
                                     allocn(nf),allocp(inferpfVL(nf, allocn(nf))), 
                                     CO2) - VLong_constraint_N(nf,allocn(nf))$NPP
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
solveVLong_root_ocn <- function(CO2) {
    fn <- function(nf) {
        photo_constraint_full_cnp(nf, inferpfVL_root_ocn(nf, allocn(nf)), 
                                  allocn(nf),allocp(inferpfVL_root_ocn(nf, allocn(nf))), 
                                  CO2) - NConsVLong_root_ocn(nf,allocn(nf))$NPP
    }
    equilnf <- uniroot(fn,interval=c(0.01,0.1))$root
    equilpf <- inferpfVL_root_ocn(equilnf, allocn(equilnf))
    equilNPP <- photo_constraint_full_cnp(equilnf, equilpf, 
                                          allocn(equilnf), allocp(equilpf), CO2)
    
    ans <- data.frame(equilnf,equilpf,equilNPP)
    return(ans)
}


# Find the long term equilibrium nf and NPP under standard conditions - by finding the root
# specifically for nuptake ~ root biomass  - O-CN approach
# i.e. N uptake as a saturating function of mineral N
solveVLong_root_gday <- function(CO2) {
    fn <- function(nf) {
        photo_constraint_full_cn(nf, allocn(nf), CO2) - NConsVLong_root_gday(nf,allocn(nf))$NPP
    }
    equilnf <- uniroot(fn,interval=c(0.001,0.1))$root
    equilNPP <- photo_constraint_full_cn(equilnf, allocn(equilnf),  CO2)
    
    ans <- data.frame(equilnf,equilNPP)
    return(ans)
}

# Find the very-long term equilibrium nf and NPP under standard conditions - by finding the root
# specifically for exudation
solveVLong_exudation <- function(CO2) {
    fn <- function(nf) {
        photo_constraint_full_cn(nf, allocn_exudation(nf), CO2) - VLong_constraint_N(nf,allocn_exudation(nf))$NPP
    }
    equilnf <- uniroot(fn,interval=c(0.01,0.1))$root
    equilNPP_N <- photo_constraint_full_cn(equilnf, 
                                            allocn_exudation(equilnf), CO2)
    
    ans <- data.frame(equilnf,equilNPP_N)
    return(ans)
}

# Find the very-long term equilibrium nf and NPP under standard conditions - by finding the root
# specifically for exudation and priming
# note leaf NC range different
solveVLong_exudation_medium <- function(CO2) {
    fn <- function(nf) {
        photo_constraint_full_cn(nf, allocn_exudation(nf), CO2) - VLong_constraint_N(nf,allocn_exudation(nf))$NPP
    }
    equilnf <- uniroot(fn,interval=c(0.004,0.1))$root
    equilNPP_N <- photo_constraint_full_cn(equilnf, 
                                           allocn_exudation(equilnf), CO2)
    
    ans <- data.frame(equilnf,equilNPP_N)
    return(ans)
}

# Find the very-long term equilibrium nf and NPP under standard conditions - by finding the root
# without priming only
solveVLong_full_cn_medium <- function(CO2) {
    fn <- function(nf) {
        photo_constraint_full_cn(nf, allocn_exudation(nf),CO2) - VLong_constraint_N(nf,allocn_exudation(nf))$NPP
    }
    equilnf <- uniroot(fn,interval=c(0.001,0.1))$root
    equilNPP <- photo_constraint_full_cn(equilnf, allocn_exudation(equilnf), CO2)
    equilpf <- "NA"
    ans <- data.frame(equilnf, "NA", equilNPP)
    return(ans)
}