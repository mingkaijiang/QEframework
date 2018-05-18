c_cost_resorb <- function(potnpp, df) {
    ### Unit of kr: kg C m-2
    
    leafC <- potnpp / sf
    
    c_cost_resorb <- kr / (leafC * df$nf)
    
    return(c_cost_resorb)
}