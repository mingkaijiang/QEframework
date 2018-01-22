
#### Functions to generate Table 1
#### Purpose:
#### To generate % diff between aCO2 and eCO2 over time under different model assumptions
#### by looping through all the simulation subfolders
#### Save a table to store all the statistics
#### Then you can make the table look nicer mannually in the manuscript
################################################################################

######################## Functions ###################################
############# print out response of CO2 into a table
doubling_effect <- function(tranDF) {
    
    # get length of tranDF
    l <- nrow(tranDF)
    
    # assign values
    outDF[i,1] <- tranDF[1,"shootn"]/tranDF[1,"shoot"]
    outDF[i,2] <- tranDF[1,"shootp"]/tranDF[1,"shoot"]
    
    outDF[i,3] <- tranDF[l,"shootn"]/tranDF[l,"shoot"]
    outDF[i,4] <- tranDF[l,"shootp"]/tranDF[l,"shoot"]
    
    outDF[i,5] <- (tranDF[6,"npp"]-tranDF[1,"npp"])/tranDF[1,"npp"] * 100.0
    outDF[i,6] <- (tranDF[16,"npp"]-tranDF[1,"npp"])/tranDF[1,"npp"] * 100.0
    outDF[i,7] <- (tranDF[106,"npp"]-tranDF[1,"npp"])/tranDF[1,"npp"] * 100.0
    outDF[i,8] <- (tranDF[l,"npp"]-tranDF[1,"npp"])/tranDF[1,"npp"] * 100.0
    
}

run_doubling_effect_table <- function() {
    #### obtain the original working directory
    cwd <- getwd()
    
    #### Setting working directory
    setwd("GDAY/analyses")
    
    #### Count number of simulations runs by counting the # folders
    dirFile <- list.dirs(path=".", full.names = TRUE, recursive = FALSE)
    
    #### Set back to the original working directory
    setwd(cwd)
    
    #### Create a dataframe to store all the numbers
    outDF <- matrix(0, ncol=8, nrow = length(dirFile))
    colnames(outDF) <- c("leafnc_350", "leafpc_350", "leafnc_700", "leafpc_700",
                         "npp1", "npp10", "npp100","npp_equil")
    outDF <- as.data.frame(outDF)
    
    #### plot continuity plot for each subdirectory
    for (i in 1:length(dirFile)) {
        ## Set file path
        FilePath <- paste(getwd(), "/GDAY/analyses/Run", i, sep="")
        
        ## Read in the eCO2 files
        inDF <- read.table(paste(FilePath, "/annual_gday_result_transient_CO2_ELE.csv", sep=""),
                           header=T,sep=",")
        
        doubling_effect(inDF)
    }
    
    write.table(outDF, paste(getwd(), "/Tables/Table1.csv", sep=""),
                col.names=T,row.names=F, sep=",")
}


######################## Program ###################################
run_doubling_effect_table()
