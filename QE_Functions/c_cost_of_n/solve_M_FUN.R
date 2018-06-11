# Find the medium term equilibrium nf and NPP under standard conditions - by finding the root
solve_M_FUN <- function(CO2,C_pass,C_slow,Nin_L) {
    fn <- function(nf) {
        photo_constraint_full(nf, alloc(nf),
                              CO2) - M_constraint_FUN(nf,alloc(nf), C_pass, C_slow, Nin_L)$NPP_grow
        
    }
    equilnf <- uniroot(fn,interval=c(0.001,0.01))$root
    
    equilNPP <- photo_constraint_full(equilnf, 
                                      alloc(equilnf), CO2)

    ans <- data.frame(equilnf,equilNPP)
    colnames(ans) <- c("nf", "NPP")
    
    return(ans)
}
