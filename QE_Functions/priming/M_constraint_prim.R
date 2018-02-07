
### Function for nutrient N constraint in medium term 
# considering priming effect on slow SOM
M_constraint_prim <- function(df, a, C_pass, C_slow, Nin_L) {
    # passed are df and a, the allocation and plant N:C ratios
    # parameters : 
    # Nin is fixed N inputs (N deposition annd fixation) in g m-2 yr-1 (could vary fixation)
    # nleach is the rate of n leaching of the mineral pool (per year)
    # Tsoil is effective soil temperature for decomposition
    # Texture is the fine soil fraction
    # ligfl and ligrl are the lignin:C fractions in the foliage and root litter
    # Cpass is the passive pool size in g C m-2
    # ncp is the NC ratio of the passive pool in g N g-1 C
    
    # set up stuffs
    len <- length(df)
    ans <- c()
    
    # burial fractions
    pas <- 0.996 - (0.85-0.68*Texture)
    psa <- 0.42
    ppa <- 0.45
    pap <- 0.004
    psp <- 0.03
    
    for (i in 1:len) {
        fPC <- function(NPP) {
            # passive and slow pool burial 
            pass <- soil_coef_prim(df[i], a[i,], NPP)
            
            # again, exclude exudation from root allocation
            omega_ap <- a[i,]$af*pass$omega_af_pass + (a[i,]$ar-a[i,]$ar*a[i,]$ariz)*pass$omega_ar_pass 
            omega_as <- a[i,]$af*pass$omega_af_slow + (a[i,]$ar-a[i,]$ar*a[i,]$ariz)*pass$omega_ar_slow 
            
            # equation for N constraint with passive, slow, wood, and leaching
            Npass <- (1-s_coef$qq_pass) * s_coef$decomp_pass * C_pass * ncp
            Nslow <- (1-s_coef$qq_slow) * s_coef$decomp_slow * C_slow * ncs
            
            U0 <- Nin_L + Npass + Nslow
            nwood <- a[i,]$aw*a[i,]$nw
            nburial <- omega_ap*ncp + omega_as*ncs
            nleach <- leachn / (1-leachn) * (a[i,]$nfl*a[i,]$af + a[i,]$nr*a[i,]$ar 
                                             + a[i,]$nw*a[i,]$aw)
            
            NPP_NC <- U0 / (nleach + nwood + nburial)
            
            # returned in kg C m-2 yr-1
            out <- NPP_NC*10^-3 - NPP 
        }
        ans[i] <- uniroot(fPC,interval=c(0.1,20), trace=T)$root
    }
    
    out <- data.frame(ans, a$ariz)
    colnames(out) <- c("NPP", "ariz")
    return(out)   
}
