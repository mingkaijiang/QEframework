#### Compute allocation coefficient of rizosphere deposition

ariz_dep_coef <- function(nf) {
    # ariz0: minimum allocation to rhizodeposition
    # ariz1: slope of allocation to rhizodeposition with leaf C:N
    # cnref: reference leaf C:N ratio for BD and NE PFTs: 25/42
    # ariz constrained to less than 0.5
    
    cnleaf <- 1.0/nf
    
    # assume you can allocate a max of 80% of root NPP as exudates
    ariz <- ariz0 + ariz1 * pmax(pmin((cnleaf - cnref)/cnref, 0.8), 0)

    # browser()
    
    return(ariz)
}
