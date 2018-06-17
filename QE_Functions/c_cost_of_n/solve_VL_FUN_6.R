# Find the very-long term equilibrium nf and NPP under standard conditions - by finding the root
solve_VL_FUN_6 <- function(CO2) {
    ### CO2 is the pre-defined CO2 concentration
    
    ### Find unique solution
    fn <- function(nf) {
        photo_constraint_FUN(nf, alloc(nf),CO2)$NPP_grow - VL_constraint_FUN_6(alloc(nf))$NPP
    }
    equil_nf <- uniroot(fn,interval=c(0.001,0.1))$root
    
    ### calculate equilibrium NPP based on equil nf ratio
    equil_NPP <- photo_constraint_FUN(equil_nf, alloc(equil_nf), CO2)$NPP_grow
    
    ans <- data.frame(equil_nf, equil_NPP)
    
    colnames(ans) <- c("nf", "NPP")
    
    return(ans)
}
