
#### Functions to generate Figure 5
#### Purpose: 
#### to draw barchart of wood, slow and passive SOM pools
#### and demonstrate the effect of wood stoichiometric flexibility
################################################################################
######### Main program
Figure_5_plotting <- function() {
    myDF <- read.csv("Tables/Stoichiometric_flexibility_table.csv")
    
    # transform the df
    temDF <- matrix(ncol=4, nrow=18)
    temDF <- as.data.frame(temDF)
    colnames(temDF) <- c("Value", "Pool", "Element", "Model")
    temDF$Pool <- rep(c("Passive", "Slow", "Wood"), each = 3)
    temDF$Element <- rep(c("C", "N", "P"), 6)
    temDF$Model <- rep(c("Variable", "Fixed"), each=9)
    temDF$Value <- c(myDF[1,2:10], myDF[2,2:10])    
    
    temDF$Pool <- as.factor(temDF$Pool)
    temDF$Element <- as.factor(temDF$Element)
    temDF$Model <- as.factor(temDF$Model)
    temDF$Value <- as.numeric(temDF$Value)
    
    ylabel <- bquote(.("Stock size g") ~ m^-2)
    
    temDF <- temDF[temDF$Element != "C", ]    
    temDF <- temDF[temDF$Element != "P", ]

    
    require(ggplot2)
    
    # making bar plots
    tiff("Plots/Figure5.tiff",
         width = 10, height = 5, units = "in", res = 300)

    p1 <- ggplot(temDF, aes(x=Element, y=Value, fill=Model)) +   
        geom_bar(position='dodge', stat='identity') +
        facet_wrap( ~ Pool, scales="free") + 
        theme(panel.background=element_blank(), axis.line = element_line(color="grey")) + 
        #theme_bw() + 
        labs(list(x = "Nutrient element", y = ylabel, fill = "Model"))
    
    print(p1)
    
    dev.off()
}



Figure_5_plotting()
