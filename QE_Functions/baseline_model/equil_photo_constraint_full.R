### NPP as function of nf and LAI (which is calculated from NPP)
### basic function: CUE dependent and based on the full model
### for cn model only
equil_photo_constraint_full <- function(nf_equ, nfdf_equ, NPP_equ, CO2_equ,
                                           rln = "walker") {
    
    ### convert NPP to g m-2 yr-1
    NPP_g_m2_yr <- NPP_equ * 1000.0
    
    ### calculate N content in canopy
    ncontent <- NPP_g_m2_yr * nfdf_equ$af / sf * nf_equ
    
    ### update sla unit from m2 kg-1 DM to m2 g-1
    sla_m2_per_g <- SLA / 1000.0
    
    ### compute LAI m2 m-2
    lai <- sla_m2_per_g * NPP_g_m2_yr * nfdf_equ$af / sf / cfrac
    
    ### Top of canopy N content
    N0 <- ncontent * kn / (1.0 - exp(-kn * lai))
    
    gamma_star <- arrh(mt, gamstar25, eag, tk)
    #gamma_star <- 32.97
    
    ### Michaelis-Menten coefficents for carboxylation by Rubisco 
    Kc <- arrh(mt, kc25, eac, tk)
    #Kc <- 234.72
    
    ### Michaelis-Menten coefficents for oxygenation by Rubisco 
    Ko <- arrh(mt, ko25, eao, tk)
    #Ko <- 216876.747
    
    ### return effective Michaelis-Menten coefficient for CO2 
    km <- (Kc * (1.0 + oi / Ko))
    #km <- 461.998
    
    ### calcualte vcmax and jmax
    if (rln == "walker") {
        log_vcmax <- 1.993 + 2.555 * log(N0) - 0.372 * log(sla_m2_per_g) + 0.422 * log(N0) * log(sla_m2_per_g)
        vcmax <- exp(log_vcmax)
        
        log_jmax <- 1.197 + 0.847 *log_vcmax
        jmax <- exp(log_jmax)
    } else if (rln == "ellsworth") {
        print("not implemented yet!")
    }
    
    ### calculate ci
    g1w <- g1 * wtfac_root
    cica <- g1w / (g1w + sqrt(vpd * PA_2_KPA))
    ci <- cica * CO2_equ
    
    ### calculate alpha: quantum efficiency
    alpha <- assim(ci, gamma_star, alpha_j/4.0, 2.0*gamma_star)
    
    ### Calculate ac, aj and asat
    ac = assim(ci, gamma_star, vcmax, km)
    aj = assim(ci, gamma_star, jmax/4.0, 2.0*gamma_star)
    asat <- pmin(aj, ac)
    
    ### This function requires further investigation
    #lue_calc <- epsilon_simplified(asat, PAR_MJ, alpha, daylen)
    lue_calc <- epsilon(asat, par, alpha, daylen)
    
    ### in umol C m-2 d-1
    lue_yr <- lue_calc * par 
    
    ### return gpp as kg m-2 yr-1
    gpp <- lue_yr * (1 - exp(-kext*SLA*nfdf_equ$af*NPP_equ/sf/cfrac)) * conv * 365 / 1000.0
    
    ##Returns G: total C production (i.e. NPP) in kg m-2 yr-1
    return( gpp * cue)
    
}