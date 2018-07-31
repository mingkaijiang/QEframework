# Find the long term equilibrium nf and NPP under standard conditions - by finding the root
solve_L_Medlyn_Dewar_no_coupling <- function(CO2,C_pass,Nin_L) {
    fn <- function(nf) {
        photo_constraint_full(nf, alloc_Medlyn_Dewar_no_coupling(nf), 
                              CO2) - L_constraint(nf,alloc_Medlyn_Dewar_no_coupling(nf),C_pass=C_pass,Nin_L=Nin_L)$NPP
    }
    equilnf <- uniroot(fn,interval=c(0.001,0.1))$root
    equilNPP <- photo_constraint_full(equilnf, 
                                    alloc_Medlyn_Dewar_no_coupling(equilnf), CO2)
    ans <- data.frame(equilnf, equilNPP)
    colnames(ans) <- c("nf", "NPP")
    
    return(ans)
}

# Find the long term equilibrium nf and NPP under standard conditions - by finding the root
solve_L_Medlyn_Dewar_linear_coupling <- function(CO2,C_pass,Nin_L) {
    fn <- function(nf) {
        photo_constraint_full(nf, alloc_Medlyn_Dewar_linear_coupling(nf), 
                              CO2) - L_constraint(nf,alloc_Medlyn_Dewar_linear_coupling(nf),C_pass=C_pass,Nin_L=Nin_L)$NPP
    }
    equilnf <- uniroot(fn,interval=c(0.001,0.1))$root
    equilNPP <- photo_constraint_full(equilnf, 
                                      alloc_Medlyn_Dewar_linear_coupling(equilnf), CO2)
    ans <- data.frame(equilnf, equilNPP)
    colnames(ans) <- c("nf", "NPP")
    
    return(ans)
}