#### setting CO2 concentrations
CO2_1 <- 400.0  # 350
CO2_2 <- 800.0  # 700

#### define parameters
nwood = 0.005     
nrho = 0.7
nretrans = 0.1
nwvar = TRUE    
kext=0.5
SLA=5
sf=0.5
cfrac = 0.45
cue = 0.5
leachn = 0.05
Nin = 0.4        
Tsoil = 15
Texture = 0.5
ligfl = 0.2
ligrl = 0.16
ncp = 0.1
PAR_MJ <- 4.0
J_2_UMOL <- 4.57
MJ_TO_J <- 1000000.0
par <- MJ_TO_J * J_2_UMOL * PAR_MJ
UMOL_TO_MOL <- 0.000001
MOL_C_TO_GRAMS_C <- 12.0
conv <- UMOL_TO_MOL * MOL_C_TO_GRAMS_C
mt <- 25.0 + 273.5  # degree to kelvin
tk <- 20.0 + 273.5  # air temperature
gamstar25 <- 42.75
eag <- 37830.0
eac <- 79430.0
eao <- 36380.0
kc25 <- 404.9
ko25 <- 278400.0
oi <- 210000.0
vpd <- 2.4
PA_2_KPA <- 0.001
wtfac_root <- 1.0
g1 <- 3.8667
alpha_j <- 0.308
daylen <- 8.0
kn <- 0.3
aroot <- 0.2
aleaf <- 0.2
ncs <- 0.01   