
#### To generate all figures
#### 
#### Note: currently the figure naming system is meaningless, 
#### The exact orders of the figures will be changed as the manuscript progresses
################################################################################

################################################################################
#### main text figures

#### Figure 1: 
#### Shoot N:C vs. production from analytical solution
#### No P considered, fixed plant stoichiometry
#### Purpose: to illustrate the concept of quasi-equilibrium interpretation
source("Plots/Figure1.R")

#### Figure 2: not plotted here, as this is the GDAY structure plot


#### Figure 3:
#### Comparing N model and NP model
source("Plots/Figure3.R")

#### Figure 4:
#### To plot shoot P:C vs. production under eCO2 and aCO2
source("Plots/Figure4.R")

#### Figure 5: 
#### Read in table 2 and draw bar chart of wood stoichiometric effect
source("Plots/Figure5.R")

#### Figure 6:
#### Explicit N uptake: barplot
#### S3: quasi-equil plot
source("Plots/Figure6_and_S3.R")

#### Figure 7:
#### comparison of GDAY vs. OCN approach
#### for N uptake 
source("Plots/Figure7.R")

#### Figure 8:
#### comparison of priming on and off
#### note: Nin changed
####       exudation turned on
source("Plots/Figure8.R")

################################################################################
#### supplementary figures

#### Figure S1:
#### To plot multi-panel transient model behaviors over the first 100 years of doubling CO2
#### comparing both fixed and variable wood stoichiometry
source("Plots/FigureS1.R")

#### Figure S2:
#### effect of wood stoichiometry
source("Plots/FigureS2.R")

