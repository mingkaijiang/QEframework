
#### Functions to generate Figure S2
#### Purpose: 
#### flexibility of wood stoichiometry
################################################################################
######### Main program
source("R/Effect_of_wood_stoichiometry.R")


wood_flex_cp_paste <- function() {
    cmd.expression <- paste0("cp Plots/Effect_of_wood_stoichiometry.tiff Plots/FigureS2.tiff")
    
    system(cmd.expression)
}


wood_flex_cp_paste()
