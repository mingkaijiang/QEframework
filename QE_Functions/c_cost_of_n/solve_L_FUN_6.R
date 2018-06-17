# Find the long term equilibrium nf and NPP under standard conditions - by finding the root
solve_L_FUN_6 <- function(CO2,C_pass,Nin_L) {
    fn <- function(nf) {
        photo_constraint_FUN(nf, alloc(nf), 
                              CO2)$NPP_grow - L_constraint_FUN_6(nf,alloc(nf),C_pass=C_pass,Nin_L=Nin_L)$NPP
    }
    equilnf <- uniroot(fn,interval=c(0.001,0.1))$root
    equilNPP <- photo_constraint_FUN(equilnf, 
                                    alloc(equilnf), CO2)$NPP_grow
    ans <- data.frame(equilnf, equilNPP)
    colnames(ans) <- c("nf", "NPP")
    
    return(ans)
}