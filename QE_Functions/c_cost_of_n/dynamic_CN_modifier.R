dynamic_CN_modifier <- function(cn_plant) {
    ### Calculate a modifier for dynamic CN ratios
    
    ### Define constants
    ### Parameters fitted to give flexible CN ranges over the perating range of N costs of the model
    acn_flex <- 0.0 # operates as the Ncost above which there is a modification in the C expenditure (to allow higher CN ratio)
    bcn_flex <- 0.0 # scalar that determines how much the C expenditure is modified for a given discrepancy between acn and the actual cost of uptake
    target_cn <- 20
    
    ### CN ratio of the plant
    delta_cn <- cn_plant - target_cn
    
    
    ### modifier
    mod <- pmax(0.0, 1.0 - (N_cost - acn_flex) / bcn_flux)
    
    
    return(mod)
}