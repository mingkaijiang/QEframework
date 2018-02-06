
# Find the medium term equilibrium nf and NPP under standard conditions - by finding the root
# specifically for exudation model, turning exudation on
solveMedium_priming <- function(CO2,Cpass,Cslow,NinL) {
    fn <- function(nf) {
        photo_constraint_full_cn(nf, allocn_exudation(nf),
                                 CO2) - NConsMedium_priming(nf,allocn_exudation(nf),Cpass,Cslow,NinL)$NPP
        
    }
    equilnf <- uniroot(fn,interval=c(0.001,0.1))$root
    equilNPP <- photo_constraint_full_cn(equilnf, 
                                         allocn_exudation(equilnf), CO2)
    ans <- data.frame(equilnf,equilNPP)
    return(ans)
}

