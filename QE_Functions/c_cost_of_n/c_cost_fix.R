c_cost_fix <- function() {
    ### units: s1: kg C kg N-1 degree C-1
    c_cost_fix <- s1 * (exp(a1+b1*Tsoil*(1-0.5*(Tsoil/c1))) - 2)
    
    return(c_cost_fix)
}