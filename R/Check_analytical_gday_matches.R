
#### Check equilibrium point matches between analytical and gday output
####
#### 
################################################################################

######################## Functions ###################################
search_long_equil_year_in_gday <- function(run, co2 = 350, equilDF, gdayDF1, gdayDF2) {
    ## Run: which run number
    ## co2 = 350 or 700
    ## equilDF is the equilibrium analytical solution df
    ## gdayDF1 is the spin up gday file
    ## gdayDF2 is the eCO2 gday file
    
    if (co2 == 350) {
        
        # obtain long term equil nc value
        equilnc <- equilDF[equilDF$Run == run & equilDF$CO2 == co2, "NC_L"]
        
        gdayDF1$shootnc <- gdayDF1$shootn/gdayDF1$shoot
        
        # search in gdayDF1 dataframe
        year1 <- which.min(abs(gdayDF1$shootnc - equilnc))
        
        # pass the year information out
        return(year1)
        
    } else if (co2 == 700) {
        
        # obtain long term equil nc value
        equilnc <- equilDF[equilDF$Run == run & equilDF$CO2 == co2, "NC_L"]
        
        shootnc <- gdayDF2$shootn/gdayDF2$shoot
        
        # search in gdayDF1 dataframe
        year2 <- which.min(abs(shootnc - equilnc))
        
        # pass the year information on of the function
        return (year2)
    }
    
}


run_check_matches <- function() {
    #### obtain the original working directory
    cwd <- getwd()
    
    #### Setting working directory
    setwd("GDAY/outputs")
    
    #### Count number of simulations runs by counting the # folders
    dirFile <- list.dirs(path=".", full.names = TRUE, recursive = FALSE)
    
    #### Set back to the original working directory
    setwd(cwd)
    
    #### create a output dataframe
    matchDF <- equilDF
    
    #### plot continuity plot for each subdirectory
    for (i in 1:length(dirFile)) {
        ## Set file path
        FilePath <- paste0("GDAY/analyses/Run", i)
        
        ## Read in the eCO2 files
        gdayDF1 <- read.table(paste0(FilePath, "/annual_gday_result_spinup.csv"),
                              header=T,sep=",")
        
        ## Read in the eCO2 files
        gdayDF2 <- read.table(paste0(FilePath, "/annual_gday_result_transient_CO2_ELE.csv"),
                              header=T,sep=",")
        print(FilePath)
        
        ## obtain file length
        l1 <- nrow(gdayDF1)
        l2 <- nrow(gdayDF2)
        
        ## calculate long-term equilibrium point
        yr1 <- search_long_equil_year_in_gday(i, co2 = 350, equilDF, gdayDF1, gdayDF2)
        yr2 <- search_long_equil_year_in_gday(i, co2 = 700, equilDF, gdayDF1, gdayDF2)
        
        
        ## add gday equilibrium values onto the dataframe
        # VL aCO2
        matchDF[matchDF$Run == i & matchDF$CO2 == 350, "NC_VL_GDAY"] <- gdayDF2[1, "shootn"]/gdayDF2[1, "shoot"]
        matchDF[matchDF$Run == i & matchDF$CO2 == 350, "PC_VL_GDAY"] <- gdayDF2[1, "shootp"]/gdayDF2[1, "shoot"]
        matchDF[matchDF$Run == i & matchDF$CO2 == 350, "NPP_VL_GDAY"] <- gdayDF2[1, "npp"]/10

        # L aCO2
        matchDF[matchDF$Run == i & matchDF$CO2 == 350, "NC_L_GDAY"] <- gdayDF1[yr1, "shootn"]/gdayDF1[yr1, "shoot"]
        matchDF[matchDF$Run == i & matchDF$CO2 == 350, "PC_L_GDAY"] <- gdayDF1[yr1, "shootp"]/gdayDF1[yr1, "shoot"]
        matchDF[matchDF$Run == i & matchDF$CO2 == 350, "NPP_L_GDAY"] <- gdayDF1[yr1, "npp"]/10
        matchDF[matchDF$Run == i & matchDF$CO2 == 350, "Equil_Year"] <- yr1
        
        
        # VL eCO2
        matchDF[matchDF$Run == i & matchDF$CO2 == 700, "NC_VL_GDAY"] <- gdayDF2[l2, "shootn"]/gdayDF2[l2, "shoot"]
        matchDF[matchDF$Run == i & matchDF$CO2 == 700, "PC_VL_GDAY"] <- gdayDF2[l2, "shootp"]/gdayDF2[l2, "shoot"]
        matchDF[matchDF$Run == i & matchDF$CO2 == 700, "NPP_VL_GDAY"] <- gdayDF2[l2, "npp"]/10

        # L eCO2
        matchDF[matchDF$Run == i & matchDF$CO2 == 700, "NC_L_GDAY"] <- gdayDF2[yr2, "shootn"]/gdayDF2[yr2, "shoot"]
        matchDF[matchDF$Run == i & matchDF$CO2 == 700, "PC_L_GDAY"] <- gdayDF2[yr2, "shootp"]/gdayDF2[yr2, "shoot"]
        matchDF[matchDF$Run == i & matchDF$CO2 == 700, "NPP_L_GDAY"] <- gdayDF2[yr2, "npp"]/10
        matchDF[matchDF$Run == i & matchDF$CO2 == 700, "Equil_Year"] <- yr2
        
        
    }
    
    write.table(matchDF, "Tables/analytical_matches_with_gday_output.csv",
                sep=",", col.names=T, row.names=F)
}

######################## Program ###################################
run_check_matches()