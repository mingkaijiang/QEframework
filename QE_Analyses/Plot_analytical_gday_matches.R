
#### Plot analytical (lines) matches with gday simulation (points)
####
#### 
################################################################################

######################## Functions ###################################
### Plot analytical solutions and gday points
plot_analytical <- function(tranDF, endyear) {
    ## library
    require(scatterplot3d)
    
    ## generate analytical solution
    Perform_Analytical_Run1(plotting = F)
    
    s3d<-scatterplot3d(out350DF$nc, out350DF$pc_VL, out350DF$NPP_350, xlim=c(0.0, 0.05),
                  ylim = c(0.0, 0.002), zlim=c(0, 3), 
                  type = "l", xlab = "Shoot N:C ratio", ylab = "Shoot P:C ratio", 
                  zlab = expression(paste("Production [kg C ", m^-2, " ", yr^-1, "]")),
                  color="cyan", lwd = 3)
    
    # constraint by nutrient and photo at CO2 = 350
    s3d$points3d(tranDF[,"shootn"]/tranDF[,"shoot"], 
                 tranDF[,"shootp"]/tranDF[,"shoot"], 
                 tranDF[,"npp"]/10, 
                 type = "p", pch = 19)
    
    # NPP constraint by very long term nutrient availability
    s3d$points3d(out350DF$nc, out350DF$pc_VL, out350DF$NPP_VL, type="l", col="tomato", lwd = 3)
    
    # equilibrated NPP for very long term nutrient and CO2 = 350
    s3d$points3d(equil350DF$nc_VL, equil350DF$pc_VL, equil350DF$NPP_VL,
                 type="h", pch = 19, col = "blue")
    
    # NPP constraint by long term nutrient availability
    s3d$points3d(out350DF$nc, out350DF$pc_VL, out350DF$NPP_350_L, type='l',col="violet", lwd = 3)
    
    
    # NPP constraint by CO2 = 700
    s3d$points3d(out700DF$nc, out700DF$pc_VL, out700DF$NPP_700, col="green", type="l", lwd = 3)
    s3d$points3d(out700DF$nc, out700DF$pc_700_L, out700DF$NPP_700, type = "l", col = "darkgreen", lwd = 3,
           lty = 3)
    
    s3d$points3d(equil350DF$nc_VL, equil350DF$pc_VL, 
                 inst700$equilNPP, type="h", col = "darkgreen", pch=19)
    
    # equilibrated NPP for very long term nutrient and CO2 = 700
    s3d$points3d(equil700DF$nc_VL, equil700DF$pc_VL, equil700DF$NPP_VL, 
                 type="h", col="orange", pch = 19)
    
    # equilibrated NPP for long term nutrient and CO2 = 700
    s3d$points3d(equil700DF$nc_L, equil700DF$pc_VL, equil700DF$NPP_L,
                 type="h", col="red", pch = 19)
    

    
    legend("topleft", c(expression(paste("Photo constraint at ", CO[2]," = 350 ppm")), 
                        expression(paste("Photo constraint at ", CO[2]," = 700 ppm")), 
                        "VL nutrient constraint", "L nutrient constraint",
                        "A", "B", "C", "D"),
           col=c("cyan","green", "tomato", "violet","blue", "darkgreen","red", "orange"), 
           lwd=c(2,2,2,2,NA,NA,NA,NA), pch=c(NA,NA,NA,NA,19,19,19,19), cex = 1.0, 
           bg = adjustcolor("grey", 0.8))
    
    
} 



### Read in gday files and analytical files and plotting
run_gday_analytical_matches <- function() {
    #### obtain the original working directory
    cwd <- getwd()
    
    #### Setting working directory
    setwd("GDAY/analyses")
    
    #### Count number of simulations runs by counting the # folders
    dirFile <- list.dirs(path=".", full.names = TRUE, recursive = FALSE)
    
    #### Set back to the original working directory
    setwd(cwd)
    
    #### set plotting folder directory
    plotPath <- "Plots/"
    
    #### plot continuity plot for each subdirectory
    for (i in 1:length(dirFile)) {
        ## Set file path
        FilePath <- paste0("GDAY/analyses/Run", i)
        
        ## Read in the eCO2 files
        gdayDF <- read.table(paste(FilePath, "/annual_gday_result_transient_CO2_ELE.csv", sep=""),
                           header=T,sep=",")
        
        print(FilePath)
        pdf(paste0("Plots/gday_analytical_matches_Run", i, ".pdf"), width=10,height=8)
        plot_analytical(gdayDF[1:4850,])
        dev.off()
    }
}

######################## Program ###################################
run_gday_analytical_matches()