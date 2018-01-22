arrh <- function(mt, k25, Ea, tk) {

        return (k25 * exp((Ea * (tk - mt)) / (mt * 8.314 * tk)));
}