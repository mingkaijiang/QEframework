# Find the medium term equilibrium nf and NPP under standard conditions - by finding the root
solve_M_FUN_6 <- function(CO2,C_pass,C_slow,Nin_L) {
    fn <- function(nf) {
        photo_constraint_FUN(nf, alloc(nf),
                              CO2)$NPP_grow - M_constraint_FUN_6(nf,alloc(nf), C_pass, C_slow, Nin_L)$NPP
        
    }
    equilnf <- uniroot(fn,interval=c(0.001,0.1))$root
    
    equilNPP <- photo_constraint_FUN(equilnf, 
                                      alloc(equilnf), CO2)$NPP_grow

    ans <- data.frame(equilnf,equilNPP)
    colnames(ans) <- c("nf", "NPP")
    
    return(ans)
}
