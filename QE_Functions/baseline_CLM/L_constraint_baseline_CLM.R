
### Function for nutrient N constraint in longterm ie passive, leaching considered
L_constraint_baseline_CLM <- function(df, a, C_pass, Nin_L) {
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
 
    decomp_pass <- 0.00013*52*soil_decomp(Tsoil) 
    
    # finding actual NPP, in g C m-2 yr-1
    for (i in 1:len) {
        fPC <- function(NPP) {
            # passive and slow pool burial 
            pass <- soil_coef_baseline_CLM(df[i], a[i,], NPP)
            
            # again, exclude exudation from root allocation
            omega_ap <- a[i,]$af*pass$omega_af_pass + a[i,]$ar*pass$omega_ar_pass 
            
            # total in N
            U0 <- Nin_L + (1-pass$qq_pass) * pass$decomp_pass * C_pass * ncp   
            
            # n burial rate for passive pool
            nburial <- omega_ap*ncp
            
            # actual NPP
            NPP_act <- (U0 - nburial * NPP) / (ncp * NPP / decomp_pass)
            
            # returned in kg C m-2 yr-1
            out <- NPP_act*10^-3 - NPP 
        }
        ans[i] <- uniroot(fPC,interval=c(0.1,20), trace=T)$root
    }
    
    out <- data.frame(ans)

    
    return(out)   
}
