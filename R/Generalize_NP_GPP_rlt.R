#### Scripts to obtain relationships between GPP and Leaf N and Leaf N, at annual timestep
####
#### For AmazonFACE and EucFACE in order to generalize results
####

#### Read in files
amazDF_N <- read.csv("~/Documents/Research/Projects/Amazon/AMAZ/drought/outputs/AmaFACE1_D_GDA_AMB_1999_2023.csv",skip=3)
amazDF_P <- read.csv("~/Documents/Research/Projects/Amazon/AMAZ/drought_p/outputs/AmaFACE1_D_GDA_AMB_1999_2023.csv",skip=3)

aDF_n <- subset(amazDF_N, YEAR == 1999)
aDF_p <- subset(amazDF_P, YEAR == 1999)

#### Create a storage DF for annual step calculation
#t.s <- seq(1999, 2023)
#aDF_n <- data.frame(t.s, NA, NA, NA, NA, NA, NA, NA, NA)
#aDF_p <- data.frame(t.s, NA, NA, NA, NA, NA, NA, NA, NA, NA)
#colnames(aDF_n) <- c("Year", "CO2", "PAR", "GPP", "CL", "NL", "LAI", "LMA", "LUE")
#colnames(aDF_p) <- c("Year", "CO2", "PAR", "GPP", "CL", "NL", "PL", "LAI", "LMA", "LUE")
#
##### Assign fluxes and stocks
#for (i in 1999:2023) {
#    aDF_n[aDF_n$Year == i, "CO2"] <- amazDF_N[amazDF_N$YEAR == i & amazDF_N$DOY == 1, "CO2"]
#    aDF_n[aDF_n$Year == i, "PAR"] <- sum(amazDF_N[amazDF_N$YEAR == i, "PAR"])
#    aDF_n[aDF_n$Year == i, "GPP"] <- sum(amazDF_N[amazDF_N$YEAR == i, "GPP"])
#    aDF_n[aDF_n$Year == i, "CL"] <- amazDF_N[amazDF_N$YEAR == i & amazDF_N$DOY == 1, "CL"]
#    aDF_n[aDF_n$Year == i, "NL"] <- amazDF_N[amazDF_N$YEAR == i & amazDF_N$DOY == 1, "NL"]
#    aDF_n[aDF_n$Year == i, "LAI"] <- amazDF_N[amazDF_N$YEAR == i & amazDF_N$DOY == 1, "LAI"]
#    aDF_n[aDF_n$Year == i, "LMA"] <- amazDF_N[amazDF_N$YEAR == i & amazDF_N$DOY == 1, "LMA"]
#    
#    aDF_p[aDF_p$Year == i, "CO2"] <- amazDF_P[amazDF_P$YEAR == i & amazDF_P$DOY == 1, "CO2"]
#    aDF_p[aDF_p$Year == i, "PAR"] <- sum(amazDF_P[amazDF_P$YEAR == i, "PAR"])
#    aDF_p[aDF_p$Year == i, "GPP"] <- sum(amazDF_P[amazDF_P$YEAR == i, "GPP"])
#    aDF_p[aDF_p$Year == i, "CL"] <- amazDF_P[amazDF_P$YEAR == i & amazDF_P$DOY == 1, "CL"]
#    aDF_p[aDF_p$Year == i, "NL"] <- amazDF_P[amazDF_P$YEAR == i & amazDF_P$DOY == 1, "NL"]
#    aDF_p[aDF_p$Year == i, "PL"] <- amazDF_P[amazDF_P$YEAR == i & amazDF_P$DOY == 1, "PL"]
#    aDF_p[aDF_p$Year == i, "LAI"] <- amazDF_P[amazDF_P$YEAR == i & amazDF_P$DOY == 1, "LAI"]
#    aDF_p[aDF_p$Year == i, "LMA"] <- amazDF_P[amazDF_P$YEAR == i & amazDF_P$DOY == 1, "LMA"]
#}
#
#### prepare variables in the correct unit
aDF_n$SLA <- 1/aDF_n$LMA
aDF_n$LEAFNC <- aDF_n$NL/aDF_n$CL

aDF_p$SLA <- 1/aDF_p$LMA
aDF_p$LEAFNC <- aDF_p$NL/aDF_p$CL
aDF_p$LEAFPC <- aDF_p$PL/aDF_p$CL

#### Compute relationships
lm.n <- lm(LUE~LEAFNC, data=aDF_n)
summary(lm.n)
# Generalized relationship for LUE (umol C umol-1 PAR) = 0.78 * LEAFNC
#                              LUE = -0.026 + 1.49 LEAFNC

lm.p <- lm(LUE~LEAFNC:LEAFPC, data=aDF_p)
summary(lm.p)
# Generalized relationship for LUE (umol C umol-1 PAR) = 1.45 * LEAFNC - 26.41 * LEAFPC
#                              LUE = -0.00099 + 467.2 * LEAFNC * LEAFPC

#### Plot to see relationships
with(aDF_n, plot(LUE~LEAFNC))

with(aDF_p, plot(LUE~LEAFNC))
with(aDF_p, plot(LUE~LEAFPC))


#### Compute relationships
lm.n <- lm(log(LUE)~-1+log(LEAFNC), data=aDF_n)
summary(lm.n)

lm.p <- lm(log(LUE)~-1+log(LEAFPC)+log(LEAFNC), data=aDF_p)
summary(lm.p)


#### Summary notes:
# There is no good relationship from simulating for AmazonFACE,
# because the range of NC and PC are very narrow
# It may make sense to just use the full photosynthesis model within the analytical solutions
# to come up with a range of nc and lue, then simplify their relationships
