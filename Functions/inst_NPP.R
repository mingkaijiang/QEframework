### To calculate instantaneous NPP to CO2 increase
inst_NPP <- function(nf, df) {
    # nf: the equil nf value at CO2 = 350 ppm
    # df: the photosynthetic constraint curve at co2 = 700 ppm
    
    # round the nf value and finding its location on df
    nf1 <- round(nf, 3)
    loc1 <- match(nf1, df[,1])
    
    # obtain its previous/subsequent locations
    if(nf1 > nf) {
        loc2 <- loc1-1
        eqDF <- rbind(df[loc2, ], df[loc1,])
    } else if (nf1 < nf) {
        loc2 <- loc1+1
        eqDF <- rbind(df[loc1, ], df[loc2,])
    } 
    
    # obtain linear equation
    lm.model <- lm(eqDF[,2]~eqDF[,1])
    
    # calculate equil NPP
    a <- as.numeric(coefficients(lm.model)[1])
    b <- as.numeric(coefficients(lm.model)[2])
    
    equilNPP <- a+b*nf
    
    if (nf1 == nf) {
        out <- df[loc1,]
    } else {
        out <- cbind(nf, equilNPP)
    }
    
    return(as.data.frame(out))
}