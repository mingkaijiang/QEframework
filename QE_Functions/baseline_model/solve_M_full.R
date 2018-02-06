

# Find the medium term equilibrium nf and NPP under standard conditions - by finding the root
# specifically for exudation model, without considering exudation
solve_M_full <- function(CO2,C_pass,C_slow,Nin_L) {
    fn <- function(nf) {
        photo_constraint_full(nf, alloc(nf),
                              CO2) - M_constraint(nf,alloc(nf), C_pass, C_slow, Nin_L)$NPP
        
    }
    equilnf <- uniroot(fn,interval=c(0.001,0.1))$root
    
    equilNPP <- photo_constraint_full(equilnf, 
                                      alloc(equilnf), CO2)

    ans <- data.frame(equilnf,equilNPP)
    colnames(ans) <- c("nf", "NPP")
    
    return(ans)
}
