#### Compute allocation coefficient of rizosphere deposition

ariz_dep_coef <- function(nf) {
    # ariz0: minimum allocation to rhizodeposition
    # ariz1: slope of allocation to rhizodeposition with leaf C:N
    # cnref: reference leaf C:N ratio for BD and NE PFTs: 25/42
    # ariz constrained to less than 0.5
    
    cnleaf <- 1.0/nf
    
    #ariz <- pmin(0.9, ariz0 + ariz1 * pmax((cnleaf - cnref)/cnref, 0))
    ariz <- ariz0 + ariz1 * pmax((cnleaf - cnref)/cnref, 0)
    
    return(ariz)
}
