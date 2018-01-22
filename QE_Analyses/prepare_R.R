
#### Prepare folder structure and libraries
####
#### 
################################################################################

#### Install libraries
if(!require(pacman))install.packages("pacman")
pacman::p_load(scatterplot3d, 
               data.table, 
               lattice, 
               plyr,
               rPython,
               grid,
               gridBase,
               ini) # add other packages needed to this list


#### Sourcing all R files in the function subdirectory
sourcefiles <- dir("QE_Functions", pattern="[.]R$", recursive = TRUE, full.names = TRUE)
for(a in sourcefiles)source(a)

#### Sourcing all QE scripts
# scriptfiles <- dir("QE_Scripts", pattern="[.]R$", recursive = TRUE, full.names = TRUE)
# for(b in scriptfiles)source(b)

#### graphic default settings
op <- par()
