CO2_fertilization_summary <- function() {
    
    inDF <- equilDF
    
    outDF <- matrix(nrow=(nrow(inDF)/2), ncol=4)
    outDF <- as.data.frame(outDF)
    names(outDF) <- c("Run", "Inst", "L", "VL")
    outDF$Run <- unique(inDF$Run)
    
    for (i in 1:nrow(outDF)) {
        outDF[i,"Inst"] <- (inDF[inDF$Run == i & inDF$CO2 == 350, "NPP_inst"] 
        - inDF[inDF$Run == i & inDF$CO2 == 350, "NPP_VL"]) / inDF[inDF$Run == i & inDF$CO2 == 350, "NPP_VL"]
        
        outDF[i,"L"] <- (inDF[inDF$Run == i & inDF$CO2 == 700, "NPP_L"] 
                            - inDF[inDF$Run == i & inDF$CO2 == 350, "NPP_VL"]) / inDF[inDF$Run == i & inDF$CO2 == 350, "NPP_VL"]
        
        outDF[i,"VL"] <- (inDF[inDF$Run == i & inDF$CO2 == 700, "NPP_VL"] 
                            - inDF[inDF$Run == i & inDF$CO2 == 350, "NPP_VL"]) / inDF[inDF$Run == i & inDF$CO2 == 350, "NPP_VL"]
    }
 
    tDF <- t(outDF[,2:4])
    colnames(tDF) <- c("Run1", "Run2", "Run3")
    tDF <- tDF * 100.0
    
    ### plot 2-d plots of nf vs. npp and nf vs. pf
    tiff("Plots/CO2_fertilization_effect.tiff",
         width = 10, height = 5, units = "in", res = 300)
    par(mar=c(5.1,6.1,2.1,2.1))
    
    barplot(tDF, xlab = "Run numbers", col=c("darkgreen", "red","orange"),
            legend=rownames(tDF), beside=T, ylim = c(0, 50), ylab="Increase in NPP [%]",
            cex.lab=1.5)
    
    
    dev.off()
}

CO2_fertilization_summary()
