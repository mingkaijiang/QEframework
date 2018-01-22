
### Calculate passive pool nc ratio
calc_passive_nc <- function (equilNPP, adf, p1 = 0.01, p2 = 1.0, nuptake) {
    # ncp: nc ratio of the passive pool
    # equilNPP: npp value
    # p1: intercept
    # p2: trend
    
    ncp <- p1 + p2 * (equilNPP * (adf$nfl*adf$af + adf$nr*(adf$ar) + adf$nw*adf$aw)) / nuptake
    
    return(ncp)
}