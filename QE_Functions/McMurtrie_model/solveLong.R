
# Find the long term equilibrium nf and NPP under standard conditions - by finding the root
solveLong <- function(co2=350,Cpass,NinL) {
    
    fn <- function(nf) {
        solvePC(nf,alloc_mcm(nf)$af,co2) - NConsLong(nf,alloc_mcm(nf),Cpass=Cpass,NinL)$NPP_N
    }
    equilnf <- uniroot(fn,interval=c(0.01,0.05))$root
    equilNPP <- solvePC(equilnf,af=alloc_mcm(equilnf)$af, co2)
    ans <- data.frame(equilnf,equilNPP)
    return(ans)
    
}