
# Find the very-long term equilibrium nf and NPP under standard conditions - by finding the root
solveVLong <- function(co2=350) {
    
    fn <- function(nf) {
        solvePC(nf,alloc_mcm(nf)$af,co2=co2) - NConsVLong(nf,alloc_mcm(nf))$NPP
    }
    equilnf <- uniroot(fn,interval=c(0.01,0.05))$root
    equilNPP <- solvePC(equilnf,af=alloc_mcm(equilnf)$af, co2=co2)
    ans <- data.frame(equilnf,equilNPP)
    return(ans)
    
}