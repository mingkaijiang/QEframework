

# Find the very-long term equilibrium nf and NPP under standard conditions - by finding the root
solve_VL_simple <- function(CO2) {
    fn <- function(nf) {
        photo_constraint_simple(nf, alloc(nf),CO2) - VL_constraint(nf,alloc(nf))$NPP
    }
    equilnf <- uniroot(fn,interval=c(0.01,0.1))$root
    equilNPP <- photo_constraint_simple(equilnf, alloc(equilnf), CO2)
    ans <- data.frame(equilnf, equilNPP)
    colnames(ans) <- c("nf", "NPP")
    return(ans)
}
