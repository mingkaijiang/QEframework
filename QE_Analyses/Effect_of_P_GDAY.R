#### Functions
P_limitation_GDAY <- function() {
    #### read in csv files
    npDF <- read.csv("GDAY/analyses/Run1/annual_gday_result_transient_CO2_ELE.csv")
    nDF <- read.csv("GDAY/analyses/Run2/annual_gday_result_transient_CO2_ELE.csv")
    
    #### Making plots
    tiff("Plots/Effect_of_P_GDAY.tiff",
         width = 8, height = 7, units = "in", res = 300)
    par(mar=c(5.1,6.1,2.1,2.1))
    
    split.screen(c(2,1))
    split.screen(c(1,2), screen=2)
    
    screen(1)
    with(nDF[1:100,], plot(npp/10~year, ylim = c(1, 3), type='l', xlab = "Year",
                           ylab = expression(paste("Production [kg C ", m^-2, " ", yr^-1, "]")),
                           col="blue", lwd = 3, cex.lab = 1.5))
    with(npDF[1:100,], points(npp/10~year, type="l", col="red", lwd = 3))
    
    legend("topright", c("N only", 
                         "N-P"),
           col=c("blue","red"), 
           lwd=c(2,2), cex = 1.2, 
           bg = adjustcolor("grey", 0.8))
    
    screen(3)
    with(nDF[1:100,], plot(nuptake*100~year, ylim = c(0, 10), type='l', xlab = "Year",
                           ylab = expression(paste("N uptake [g N ", m^-2, " ", yr^-1, "]")),
                           col="blue", lwd = 3, cex.lab = 1.0))
    with(npDF[1:100,], points(nuptake*100~year, type="l", col="red", lwd = 3))
    
    screen(4)
    with(npDF[1:100,], plot(puptake*100~year, ylim = c(0.2, 0.4), type='l', xlab = "Year",
                           ylab = expression(paste("P uptake [g P ", m^-2, " ", yr^-1, "]")),
                           col="red", lwd = 3, cex.lab = 1.0))
    
    close.screen(all=T)
    
    dev.off()
    
}


#### Program
P_limitation_GDAY()
