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

asat <- 0.29
PAR_MJ <- 4.0
J_2_UMOL <- 4.57
MJ_TO_J <- 1000000.0
par <- MJ_TO_J * J_2_UMOL * PAR_MJ
alpha <- 0.06047569
daylen <- 8
kext <- 0.5

asat_range <- seq(0.001, 50, by = 0.001)
lue_out1 <- epsilon(asat_range, par, alpha, daylen)

PAR_MJ <- seq(1, 16, by = 1)
par_range <- MJ_TO_J * J_2_UMOL * PAR_MJ

lue_out2 <- epsilon(asat, par_range, alpha, daylen)

outDF1 <- data.frame(asat_range, lue_out1)
outDF2 <- data.frame(PAR_MJ, lue_out2)
outDF1$PAR_MJ <- 4.0
outDF1$alpha <- alpha
outDF2$asat_range <- 0.29
outDF2$alpha <- alpha


alpha_range <- seq(0.01, 0.2, by = 0.01)
lue_out3 <- epsilon(asat, par, alpha_range, daylen)

daylen_range <- seq(2, 14, by = 1)
lue_out4 <- epsilon(asat, par, alpha, daylen_range)

outDF1 <- data.frame(asat_range, 4.0, alpha, daylen, lue_out1)
outDF2 <- data.frame(asat, PAR_MJ, alpha, daylen, lue_out2)
outDF3 <- data.frame(asat, 4.0, alpha_range, daylen, lue_out3)
outDF4 <- data.frame(asat, 4.0, alpha, daylen_range, lue_out4)
colnames(outDF1) <- c("asat", "par", "alpha", "daylen", "lue")
colnames(outDF2) <- c("asat", "par", "alpha", "daylen", "lue")
colnames(outDF3) <- c("asat", "par", "alpha", "daylen", "lue")
colnames(outDF4) <- c("asat", "par", "alpha", "daylen", "lue")


outDF <- rbind(outDF1, outDF2, outDF3, outDF4)

fit.equ <- lm(lue ~ asat + par + alpha + daylen, data=outDF)
summary(fit.equ)




