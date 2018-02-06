# Find the medium term equilibrium nf and NPP under standard conditions - by finding the root
# for priming
solve_M_full_prim <- function(CO2,C_pass,C_slow,Nin_L) {
    fn <- function(nf) {
        photo_constraint_full(nf, alloc_prim(nf),
                              CO2) - M_constraint_prim(nf,alloc_prim(nf), C_pass, C_slow, Nin_L)$NPP
        
    }
    equilnf <- uniroot(fn,interval=c(0.001,0.1))$root
    
    equilNPP <- photo_constraint_full(equilnf, 
                                      alloc_prim(equilnf), CO2)

    ans <- data.frame(equilnf,equilNPP)
    colnames(ans) <- c("nf", "NPP")
    
    return(ans)
}
