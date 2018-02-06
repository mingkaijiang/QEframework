# Find the long term equilibrium nf and NPP under standard conditions - by finding the root
solve_L_full_root_gday <- function(CO2,C_pass,Nin_L) {
    fn <- function(nf) {
        photo_constraint_full(nf, alloc(nf), 
                              CO2) - L_constraint_root_gday(nf,alloc(nf),
                                                            C_pass=C_pass,Nin_L=Nin_L)$NPP
    }
    equilnf <- uniroot(fn,interval=c(0.001,0.1))$root
    equilNPP <- photo_constraint_full(equilnf, 
                                    alloc(equilnf), CO2)
    ans <- data.frame(equilnf, equilNPP)
    colnames(ans) <- c("nf", "NPP")
    
    return(ans)
}