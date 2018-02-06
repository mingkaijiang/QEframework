# Find the very-long term equilibrium nf and NPP under standard conditions - by finding the root
solve_VL_full_prim <- function(CO2) {
    ### CO2 is the pre-defined CO2 concentration
    
    ### Find unique solution
    fn <- function(nf) {
        photo_constraint_full(nf, alloc_prim(nf),CO2) - VL_constraint(nf,alloc_prim(nf))$NPP
    }
    equil_nf <- uniroot(fn,interval=c(0.001,0.1))$root
    
    ### calculate equilibrium NPP based on equil nf ratio
    equil_NPP <- photo_constraint_full(equil_nf, alloc_prim(equil_nf), CO2)
    
    ans <- data.frame(equil_nf, equil_NPP)
    
    colnames(ans) <- c("nf", "NPP")
    
    return(ans)
}
