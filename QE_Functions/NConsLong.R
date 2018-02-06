

### Function for nutrient N constraint in longterm ie passive, leaching, wood considered
# specifically for exudation
NConsLong_exudation <- function(df, a, Cpass, NinL) {
    # passed are df and a, the allocation and plant N:C ratios
    # parameters : 
    # Nin is fixed N inputs (N deposition annd fixation) in g m-2 yr-1 (could vary fixation)
    # nleach is the rate of n leaching of the mineral pool (per year)
    # Tsoil is effective soil temperature for decomposition
    # Texture is the fine soil fraction
    # ligfl and ligrl are the lignin:C fractions in the foliage and root litter
    # Cpass is the passive pool size in g C m-2
    # ncp is the NC ratio of the passive pool in g N g-1 C
    
    len <- length(df)
    
    ans <- c()
    
    for (i in 1:len) {
        fPC <- function(NPP) {
            # passive pool burial 
            pass <- passive_exudation(df[i], a[i,], NPP)   
            omegap <- a[i,]$af*pass$omegaf + a[i,]$ar*pass$omegar 
            
            # equation for N constraint with passive, wood, and leaching
            U0 <- NinL + (1-pass$qq) * pass$decomp * Cpass * ncp   
            nwood <- 0 # a[i,]$aw*a[i,]$nw
            nburial <- omegap*ncp
            nleach <- leachn / (1-leachn) * (a[i,]$nfl*a[i,]$af + a[i,]$nr*a[i,]$ar 
                                             + a[i,]$nw*a[i,]$aw)
            
            NPP_NC <- U0 / (nleach + nwood + nburial)
            
            out <- NPP_NC*10^-3 - NPP # returned in kg C m-2 yr-1
        }
        
        ans[i] <- uniroot(fPC,interval=c(0.1,20), trace=T)$root
        
    }

    out <- data.frame(ans, a$aw)
    colnames(out) <- c("NPP", "aw")
    return(out)   
}

### Function for nutrient N constraint in longterm ie passive, leaching, wood considered
# specifically for exudation
NConsLong_exudation_medium <- function(df, a, Cpass, NinL) {
    # passed are df and a, the allocation and plant N:C ratios
    # parameters : 
    # Nin is fixed N inputs (N deposition annd fixation) in g m-2 yr-1 (could vary fixation)
    # nleach is the rate of n leaching of the mineral pool (per year)
    # Tsoil is effective soil temperature for decomposition
    # Texture is the fine soil fraction
    # ligfl and ligrl are the lignin:C fractions in the foliage and root litter
    # Cpass is the passive pool size in g C m-2
    # ncp is the NC ratio of the passive pool in g N g-1 C
    
    len <- length(df)
    
    ans <- c()
    
    for (i in 1:len) {
        fPC <- function(NPP) {
            # passive pool burial 
            pass <- slow_exudation(df[i], a[i,], NPP)   
            omegap <- a[i,]$af*pass$omegafp + a[i,]$ar*pass$omegarp 
            
            # equation for N constraint with passive, wood, and leaching
            U0 <- NinL + (1-pass$qpq) * pass$decomp_p * Cpass * ncp   
            nwood <- 0 # a[i,]$aw*a[i,]$nw
            nburial <- omegap*ncp
            nleach <- leachn / (1-leachn) * (a[i,]$nfl*a[i,]$af + a[i,]$nr*a[i,]$ar 
                                             + a[i,]$nw*a[i,]$aw)
            
            NPP_NC <- U0 / (nleach + nwood + nburial)
            
            out <- NPP_NC*10^-3 - NPP # returned in kg C m-2 yr-1
        }
        
        ans[i] <- uniroot(fPC,interval=c(0.1,20), trace=T)$root
        
    }
    
    out <- data.frame(ans, a$aw)
    colnames(out) <- c("NPP", "aw")
    return(out)   
}
