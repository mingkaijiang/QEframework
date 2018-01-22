#### Function to adjust residence time of the passive pool
adjust_passive_residence_time <- function(df, a, active_in) {
    
    # compute active out at daily timestep
    # c_into_exud <- npp * a$ar * a$ariz
    
    # try to use c_from_exud_into_active as a proxy for c out from active pool
    # c_from_exud_into_active <- c_into_exud * cue_mic
    
    # residence time, suppose active_in is active_out
    rt_pass_pool = (1.0 / prime_y) / 
    pmax(0.01, (active_in / (active_in + prime_z)))
    
    # compute decomposition rate
    kdec7_new = 1.0 / rt_pass_pool;
    
    # annual timestep
    kdec7 <- kdec7_new 
    
    return(kdec7)
}

#### Function to adjust residence time of the slow pool
#adjust_slow_residence_time <- function(df, a, active_in) {
#    
#    # residence time, suppose active_in is active_out
#    rt_slow_pool = (1.0 / prime_y_slow) / 
#        pmax(0.01, (active_in / (active_in + prime_z_slow)))
#    
#    # compute decomposition rate
#    kdec6_new = 1.0 / rt_slow_pool;
#    
#    # annual timestep
#    kdec6 <- kdec6_new 
#    
#    return(kdec6)
#}#