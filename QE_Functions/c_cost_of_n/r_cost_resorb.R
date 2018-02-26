c_cost_resorb <- function(df) {
    ### Unit of kr: kg C m-2
    
    c_cost_resorb <- kr / df$nf
    
    return(c_cost_resorb)
}