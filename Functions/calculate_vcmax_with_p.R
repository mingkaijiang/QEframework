calculate_vcmax_with_p <- function(tk, N0, P0, mt) {
    log_vcmax = 3.946 + 0.921 * log(N0) + 0.121 * log(P0) + 0.282 * log(N0) * log(P0)
    vcmax = exp(log_vcmax)
    
    return(vcmax)
}