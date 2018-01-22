calculate_quantum_efficiency <- function(ci, gamma_star) {
    return(assim(ci, gamma_star, alpha/4.0, 2.0*gamma_star))
}