
# Find the medium term equilibrium nf and NPP under standard conditions - by finding the root
# specifically for original model (i.e. Figure 1)
solveMedium <- function(CO2,Cpass,Cslow,NinL) {
    fn <- function(nf) {
        solvePC(nf,alloc_mcm(nf)$af,CO2) - NConsMedium(nf,alloc_mcm(nf),Cpass,Cslow,NinL)$NPP
        
    }
    equilnf <- uniroot(fn,interval=c(0.01,0.05))$root
    equilNPP <- solvePC(equilnf,af=alloc_mcm(equilnf)$af, CO2)
    
    ans <- data.frame(equilnf,equilNPP)
    return(ans)
}

