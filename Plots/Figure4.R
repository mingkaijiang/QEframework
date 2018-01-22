
#### Functions to generate Figure 4
#### Purpose: 
#### To demonstrate how shoot PC changes with NC
#### i.e. simply repeat analytical run 1 plot

################################################################################
######### Main program

analytical_1_cp_paste <- function() {
    cmd.expression <- paste0("cp Plots/Analytical_Run1_2d.tiff Plots/Figure4.tiff")
        
    system(cmd.expression)
}


analytical_1_cp_paste()
