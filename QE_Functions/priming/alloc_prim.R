
alloc_prim <- function(nf) {
    ### Allocation and plant N concentrations - required for both PS constraint and NC constraint
    ### specifically for exudation and priming
    # parameters
    # nf is the NC ratio of foliage
    # nw is the NC ratio of wood if fixed; otherwise the ratio of wood N:C to foliage N:C
    # nwvar is whether or not to allow wood NC to vary
    # nrho is the ratio of root N:C to foliage N:C
    # nretrans is the fraction of foliage N:C retranslocated  

    len <- length(nf)
    ar <- af <- aw <- ariz <- nw <- nr <- nriz <- rep(0,len)  # initialise
    for (i in 1:len) {
        ar[i] <- aroot
        af[i] <- aleaf
        aw[i] <- 1 - ar[i] - af[i]
    }
    
    # N concentrations of rest of plant   # in g N g-1 C
    if (nwvar == FALSE) {
        nw <- nwood
    } else {
        nw <- nwood*nf 
    }
    
    nfl <- (1.0-nretrans)*nf     
    nr <- nrho*nf
    nriz <- nr

    # update allocation coefficient for rhizodeposition
    ariz <- ariz_dep_coef(nf)
    
    ret <- data.frame(nf,nfl,nw,nr,nriz,af,aw,ar,ariz)
    return(ret)
}
