### LUE function of N & Ca
LUE <- function(nf, co2, LUE0, Nref) {
    
    CaResp <- 1.632 * (co2-60.9) / (co2+121.8)
    Nresp <- min(nf/Nref, 1)
    return(LUE0 * CaResp * Nresp)
}
