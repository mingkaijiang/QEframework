epsilon <- function(asat, par, alpha, daylen) {
    # Canopy scale LUE using method from Sands 1995, 1996;
    # Parameters
    # ----------
    # asat: light-saturated photosynthetic rate at the top of the canopy
    # par: photosynthetically active radiation (umol m-2 d-1)
    # theta: curvature of photosynthetic light response curve
    # alpha: quantum yield of photosynthesis (mol mol-1)
    # Returns:
    # --------
    # LUE: integrated light use efficiency over the canopy (umol C umol -1 PAR)
    
    # subintervals scalar, i.e. 6 intervals 
    delta <- 0.16666666667
    #delta <- 0.0833333
        
    # number of seconds of daylight 
    h <- daylen * 3600.0;
    
    # theta
    theta <- 0.7
            
    # normalised daily irradiance 
    q <- pi * kext * alpha * par / (2 * h * asat);
    integral_g <- 0.0
    
    for (i in c(1,3,5,7,9,11)) {
    sinx <- sin(pi * i / 12.)
    arg1 <- sinx
    arg2 <- 1.0 + q * sinx
    arg3 <- (sqrt(((1.0 + q * sinx)^2) - 4.0 * theta * q * sinx))
    integral_g <- integral_g + (arg1 / (arg2 + arg3))
    }
    integral_g = delta * integral_g
    lue = alpha * integral_g * pi

    return (lue)
}

epsilon_simplified <- function(asat, par, alpha, daylen) {
    # simplified function for computing LUE, for using annual results
    lue <- 2.595e-2 + 8.955e-4 * asat - 1.53e-3 * par - 1.118e-1 * alpha
    return(lue)
    
}