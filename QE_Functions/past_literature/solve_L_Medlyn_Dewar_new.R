# Find the long term equilibrium nf and NPP under standard conditions - by finding the root
solve_L_Medlyn_Dewar_no_coupling_new <- function(CO2) {
    fn <- function(nf) {
        photo_constraint_full(nf, alloc_Medlyn_Dewar_no_coupling(nf), 
                              CO2) - L_constraint_Medlyn_Dewar_no_coupling(alloc_Medlyn_Dewar_no_coupling(nf))$NPP
    }
    equilnf <- uniroot(fn,interval=c(0.001,0.1))$root
    equilNPP <- photo_constraint_full(equilnf, 
                                    alloc_Medlyn_Dewar_no_coupling(equilnf), CO2)
    ans <- data.frame(equilnf, equilNPP)
    colnames(ans) <- c("nf", "NPP")
    
    return(ans)
}

# Find the long term equilibrium nf and NPP under standard conditions - by finding the root
solve_L_Medlyn_Dewar_linear_coupling_new <- function(CO2) {
    fn <- function(nf) {
        photo_constraint_full(nf, alloc_Medlyn_Dewar_linear_coupling(nf), 
                              CO2) - L_constraint_Medlyn_Dewar_linear_coupling(alloc_Medlyn_Dewar_linear_coupling(nf))$NPP
    }
    equilnf <- uniroot(fn,interval=c(0.001,0.1))$root
    equilNPP <- photo_constraint_full(equilnf, 
                                      alloc_Medlyn_Dewar_linear_coupling(equilnf), CO2)
    ans <- data.frame(equilnf, equilNPP)
    colnames(ans) <- c("nf", "NPP")
    
    return(ans)
}