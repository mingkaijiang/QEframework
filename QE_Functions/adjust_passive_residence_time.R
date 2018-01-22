#### Function to adjust residence time of the passive pool
adjust_passive_residence_time <- function(df, a, active_in) {
    
    # residence time, suppose active_in is active_out
    rt_pass_pool = (1.0 / prime_y) / 
    pmax(0.01, (active_in / (active_in + prime_z)))
    
    # compute decomposition rate
    kdec7_new = 1.0 / rt_pass_pool;
    
    # annual timestep
    kdec7 <- kdec7_new 
    
    return(kdec7)
}
