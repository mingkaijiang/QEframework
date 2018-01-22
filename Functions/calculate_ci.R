calculate_ci <- function(vpd, Ca) {
    g1w <- g1 * wtfac_root
    cica = g1w / (g1w + sqrt(vpd * PA_2_KPA))
    ci = cica * Ca
    return(ci)
}