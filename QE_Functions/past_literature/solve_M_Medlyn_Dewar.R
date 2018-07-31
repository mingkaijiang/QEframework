# Find the medium term equilibrium nf and NPP under standard conditions - by finding the root
solve_M_Medlyn_Dewar_no_coupling <- function(CO2,C_pass,C_slow,Nin_L) {
    fn <- function(nf) {
        photo_constraint_full(nf, alloc_Medlyn_Dewar_no_coupling(nf),
                              CO2) - M_constraint(nf,alloc_Medlyn_Dewar_no_coupling(nf), C_pass, C_slow, Nin_L)$NPP
        
    }
    equilnf <- uniroot(fn,interval=c(0.001,0.1))$root
    
    equilNPP <- photo_constraint_full(equilnf, 
                                      alloc_Medlyn_Dewar_no_coupling(equilnf), CO2)

    ans <- data.frame(equilnf,equilNPP)
    colnames(ans) <- c("nf", "NPP")
    
    return(ans)
}

solve_M_Medlyn_Dewar_linear_coupling <- function(CO2,C_pass,C_slow,Nin_L) {
    fn <- function(nf) {
        photo_constraint_full(nf, alloc_Medlyn_Dewar_linear_coupling(nf),
                              CO2) - M_constraint(nf,alloc_Medlyn_Dewar_linear_coupling(nf), C_pass, C_slow, Nin_L)$NPP
        
    }
    equilnf <- uniroot(fn,interval=c(0.001,0.1))$root
    
    equilNPP <- photo_constraint_full(equilnf, 
                                      alloc_Medlyn_Dewar_linear_coupling(equilnf), CO2)
    
    ans <- data.frame(equilnf,equilNPP)
    colnames(ans) <- c("nf", "NPP")
    
    return(ans)
}
