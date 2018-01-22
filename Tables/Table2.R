
#### Functions to generate Table on priming effect summary
################################################################################

######################## Functions ###################################
compute_Table2 <- function(destDir) {

   
    
    ######## Save the output table
    write.csv(out.tab, paste0(destDir, "/Table2.csv"), row.names=F)
    
}


######################## Script ###################################
compute_Table2(destDir="Tables")