#### Check continuity in transient files
#### For both pools and fluxes
#### Can set end year
#### Different from Check_continuity in that it doesn't start from spin-up files
#### And it includes fluxes in addition to pools
################################################################################


######################## Main program ###################################
## Read in spin up, aCO2 and eCO2 files
npDF <- read.table( "GDAY/analyses/Run1/annual_gday_result_transient_CO2_ELE.csv",
                   header=T,sep=",")
nDF <- read.table( "GDAY/analyses/Run2/annual_gday_result_transient_CO2_ELE.csv",
                    header=T,sep=",")

## Plot N vs. NP comparison of NPP, Nuptake and Puptake
tiff("Plots/FigureS1.tiff",
     width = 10, height = 4, units = "in", res = 300)
par(mfrow=c(1,3), mar=c(5.1,5.1,2.1,2.1))

# NPP
with(npDF[1:100,], plot(npp*100~year, type="l", 
                        ylab = expression(paste("NPP [g C ", m^-2, " ", yr^-1, "]", sep="")),
                        xlab = "Year",  lwd = 2.0, 
                        ylim = c(1500, 1900), cex.lab = 1.5))
with(nDF[1:100,], points(npp*100~year, type="l", lty = 2, lwd = 2.0))

arrows(5, 1500, 5, 1520, length=0.1)

# Nuptake
with(npDF[1:100,], plot(nuptake*100~year, type="l", 
                        ylab = expression(paste("N uptake [g N ", m^-2, " ", yr^-1, "]", sep="")),
                        xlab = "Year",  lwd = 2.0, ylim = c(0, 8), cex.lab = 1.5))
with(nDF[1:100,], points(nuptake*100~year, type="l", lty = 2, lwd = 2.0))
arrows(5, 8, 5, 7.6, length=0.1)

# Puptake
with(npDF[1:100,], plot(puptake*100~year, type="l", 
                        ylab = expression(paste("P uptake [g P ", m^-2, " ", yr^-1, "]", sep="")),
                        xlab = "Year",  lwd = 2.0, ylim = c(0.2, 0.4), cex.lab = 1.5))
arrows(5, 0.35, 5, 0.33, length=0.1)


legend("topright", c("N-P model", "N model"),
       lty = c(1, 2), lwd = 1.2, cex = 1.5)

dev.off()

rm(nDF, npDF)

### subset a period only
#tranDF <- inDF[1:100,]
#
### add extra variables
#tranDF$shootcn <- tranDF$shoot / tranDF$shootn
#tranDF$shootcp <- tranDF$shoot / tranDF$shootp
#tranDF$stemcn <- tranDF$stem / tranDF$stemn
#tranDF$stemcp <- tranDF$stem / tranDF$stemp
#tranDF$rootcn <- tranDF$root / tranDF$rootn
#tranDF$rootcp <- tranDF$root / tranDF$rootp
#
### Plotting
#tiff("Plots/Figure_transient_variable_wood_stoichiometry.tiff",
#     width = 10, height = 7, units = "in", res = 300)
#par(mar=c(5.1,5.1,2.1,2.1))
#
##### Fluxes
### Create multi-panel dataset
#panelDF <- melt(subset(tranDF, select=c("year", "npp", "shoot", "shootcn", "shootcp",
#                                        "nuptake", "stem", "stemcn", "stemcp",
#                                        "puptake", "root", "rootcn", "rootcp")),
#                id.var="year")
#
### Plot
#xyplot(value~year|variable, data=panelDF, type="l", lwd = 2.5, col = "blue",
#             scales=list(y=list(relation="free")),
#             layout=c(4,3), ylab = "", xlab = "Year")
#
#
#dev.off()
#
#########################
#
### Read in spin up, aCO2 and eCO2 files
#inDF <- read.table( "GDAY/analyses/Run3/annual_gday_result_transient_CO2_ELE.csv",
#                    header=T,sep=",")
#
### subset a period only
#tranDF <- inDF[1:100,]
#
### add extra variables
#tranDF$shootcn <- tranDF$shoot / tranDF$shootn
#tranDF$shootcp <- tranDF$shoot / tranDF$shootp
#tranDF$stemcn <- tranDF$stem / tranDF$stemn
#tranDF$stemcp <- tranDF$stem / tranDF$stemp
#tranDF$rootcn <- tranDF$root / tranDF$rootn
#tranDF$rootcp <- tranDF$root / tranDF$rootp
#
#
### Plotting
#tiff("Plots/Figure_transient_fixed_wood_stoichiometry.tiff",
#     width = 10, height = 7, units = "in", res = 300)
#par(mar=c(5.1,5.1,2.1,2.1))
#
##### Fluxes
### Create multi-panel dataset
#panelDF <- melt(subset(tranDF, select=c("year", "npp", "shoot", "shootcn", "shootcp",
#                                        "nuptake", "stem", "stemcn", "stemcp",
#                                        "puptake", "root", "rootcn", "rootcp")),
#                id.var="year")
#
### Plot
#xyplot(value~year|variable, data=panelDF, type="l", lwd = 2.5, col = "blue",
#       scales=list(y=list(relation="free")),
#       layout=c(4,3), ylab = "", xlab = "Year",
#       ylim = list(c(9, 12), c(3.7, 4.2), c(40, 70), c(1390, 1800),
#                c(0.06, 0.08), c(572, 575), c(190, 210), c(3320, 3340), 
#                c(0.0029, 0.0033), c(1.2, 1.5), c(65, 100), c(2000, 2900)))
#
#
#dev.off()
