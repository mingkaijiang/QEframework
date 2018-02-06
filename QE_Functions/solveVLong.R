
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