dynamic_CN_modifier <- function(cn_plant, N_cost) {
    ### Calculate a modifier for dynamic CN ratios
    
    ### Define constants
    ### Parameters fitted to give flexible CN ranges over the perating range of N costs of the model
    acn_flex <- 0.0 # operates as the Ncost above which there is a modification in the C expenditure (to allow higher CN ratio)
    bcn_flex <- 1 # scalar that determines how much the C expenditure is modified for a given discrepancy between acn and the actual cost of uptake
    ccn_flex <- 150.0
    target_cn <- 20
    
    ### modifier - response of C expenditure to N uptake cost
    mod <- pmax(0.0, 1.0 - (N_cost - acn_flex) / bcn_flex)
    
    ### CN ratio of the plant
    delta_cn <- cn_plant - target_cn
    
    mod2 <- c()
    
    ### modifier - response of C expenditure to plant C:N ratios
    for (k in 1:length(delta_cn)) {
        if (delta_cn[k] > 0) {
            mod2[k] <- mod[k] + 0.5 * (delta_cn[k]/ccn_flex)
        } else if (delta_cn[k] < 0) {
            mod2[k] <- mod[k] + (1 - mod[k]) * pmin(1, delta_cn[k]/ccn_flex)
        }
    }

    ### prevent unrealistically high CN ratios
    mod3 <- pmax(pmin(1.0, mod2), 0.5)
    
    return(mod3)
}