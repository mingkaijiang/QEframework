
### Allocation and plant N concentrations - required for both PS constraint and NC constraint
alloc_mcm <- function(nf,  nwvar = FALSE, nwood = 0.005, rho = 0.7, retrans = 0.5) {
    # parameters
    # nf is the NC ratio of foliage
    # nw is the NC ratio of wood if fixed; otherwise the ratio of wood N:C to foliage N:C
    # nwvar is whether or not to allow wood NC to vary
    # rho is the ratio of root N:C to foliage N:C
    # retrans is the fraction of foliage N:C retranslocated  
    
    len <- length(nf)
    ar <- af <- aw <- nw <- nr <- rep(0,len)  # initialise
    for (i in 1:len) {
        ar[i] <- 0.2
        af[i] <- 0.2
        aw[i] <- 1 - ar[i] - af[i]
    }
    
    # N concentrations of rest of plant   # in g N g-1 C
    if (nwvar == FALSE) {
        nw <- nwood
    } else {
        nw <- nwood*nf 
    }
    rho <- 0.7
    nr <- rho*nf
    nfl <- (1-retrans)*nf
    ret <- data.frame(nf,nfl,nw,nr,af,aw,ar)
    return(ret)
}
