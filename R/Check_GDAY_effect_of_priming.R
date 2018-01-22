#### Check GDAY simulated effect of root exudation on nmineralization and related production

myDF1 <- read.csv("GDAY/outputs/Run1/Quasi_equil_transient_CO2_AMB.csv",skip=1)
myDF2 <- read.csv("GDAY/outputs/Run9/Quasi_equil_transient_CO2_AMB.csv",skip=1)

# fluxes
with(myDF1, plot(nmineralisation))
with(myDF2, points(nmineralisation, col="red"))

with(myDF1, plot(nuptake, ylim=c(0.000204, 0.000206)))
with(myDF2, points(nuptake, col="red"))

with(myDF1, plot(puptake, ylim=c(0.0000086, 0.0000088)))
with(myDF2, points(puptake, col="red"))

plot(myDF1$nuptake, myDF2$nuptake, xlab = "exudation off",
     ylab = "exudation on")
abline(a=0,b=1, col="red")

with(myDF1, plot(npp))
with(myDF2, points(npp, col="red"))

# stocks
with(myDF1, plot(activesoil, ylim=c(0, 5)))
with(myDF2, points(activesoil, col="red"))

with(myDF1, plot(slowsoil, ylim=c(0, 50)))
with(myDF2, points(slowsoil, col="red"))

with(myDF1, plot(passivesoil, ylim=c(20, 50)))
with(myDF2, points(passivesoil, col="red"))


with(myDF1, plot(activesoiln, ylim=c(0, 0.2)))
with(myDF2, points(activesoiln, col="red"))

with(myDF1, plot(activesoilp, ylim=c(0, 0.01)))
with(myDF2, points(activesoilp, col="red"))

with(myDF1, plot(slowsoiln, ylim=c(0, 5)))
with(myDF2, points(slowsoiln, col="red"))

with(myDF1, plot(passivesoiln, ylim=c(0, 5)))
with(myDF2, points(passivesoiln, col="red"))

with(myDF1, plot(inorgn, ylim=c(0.078, 0.079)))
with(myDF2, points(inorgn, col="red"))


# CN ratios
with(myDF1, plot(activesoil/activesoiln, ylim=c(10, 20)))
with(myDF2, points(activesoil/activesoiln, col="red"))

with(myDF1, plot(slowsoil/slowsoiln, ylim=c(10, 30)))
with(myDF2, points(slowsoil/slowsoiln, col="red"))

with(myDF1, plot(passivesoil/passivesoiln, ylim=c(0, 15)))
with(myDF2, points(passivesoil/passivesoiln, col="red"))


# plant allocation
with(myDF1, plot(root/activesoiln, ylim=c(10, 20)))
with(myDF2, points(activesoil/activesoiln, col="red"))

## Comments:
## With Exudation turned on, we have higher NPP, higher N uptake,
##                           lower active, slow and passive soil C content
##                           lower active, slow and passive soil N content
## Root exudated C enters active SOM;
## To meet the CN demand of the active SOM pool, we need extra N from N mineralization;
## More N mineralized, bigger N mineral pool size, and higher N uptake;
## higher N uptake, higher NPP;
## Faster turnaround of the slow pool, less time to accumulate C, less slow pool;
## Less slow pool size, less active and passive pool


#### turn priming off in gday simulation
myDF2 <- read.csv("GDAY/outputs/Run9/Quasi_equil_transient_CO2_AMB.csv",skip=1)

tsoil <- 15
kdec5 <- 7.305 / 365.0
kdec6 <- 0.198279 / 365.0
kdec7 <- 0.006783 / 365.0
finesoil <- 0.5
soiltext <- 1.0 - (0.75 - finesoil)
tempfunc <- max(0.0, 0.0326 + 0.00351 * (tsoil^1.652) - ((tsoil / 41.748)^7.19))
d5prime <- kdec5 * (1 - 0.75 * soiltext) * tempfunc
d6prime <- kdec6 * tempfunc
d7prime <- kdec7 * tempfunc
    
# find pas coefficient
pas <- myDF2$active_to_slow / d5prime / myDF2$activesoil
soiltext2 <- (0.85 -  (0.996 - pas)) / 0.68
finesoil2 <- 0.75 - (1 - soiltext2)

# find psa coefficient; original value = 0.42
psa <- myDF2$slow_to_active / d6prime / myDF2$slowsoil

# find pap coefficient; original value = 0.004
pap <- myDF2$active_to_passive / d5prime / myDF2$activesoil

# find psp coefficient; original value = 0.03
psp <- myDF2$slow_to_passive / d6prime / myDF2$slowsoil

# find ppa coefficient; original value = 0.45
ppa <- myDF2$passive_to_active / d7prime / myDF2$passivesoil

# save a df
outDF <- data.frame(c("pas", "psa", "pap", "ppa", "psp", "kdec6"), NA, NA, NA)
colnames(outDF) <- c("coef", "analytical", "prime_off", "prime_on")
outDF[1,"analytical"] <- 0.996 - (0.85 - 0.68 * soiltext)
outDF[2,"analytical"] <- 0.42
outDF[3,"analytical"] <- 0.004
outDF[4,"analytical"] <- 0.45
outDF[5,"analytical"] <- 0.03
outDF[6,"analytical"] <- kdec6

outDF[1,"prime_off"] <- pas[1]
outDF[2,"prime_off"] <- psa[1]
outDF[3,"prime_off"] <- pap[1]
outDF[4,"prime_off"] <- ppa[1]
outDF[5,"prime_off"] <- psp[1]
outDF[6,"prime_off"] <- kdec6

### turn priming on and check results
myDF2 <- read.csv("GDAY/outputs/Run9/Quasi_equil_transient_CO2_AMB.csv",skip=1)

# find pas coefficient
pas <- myDF2$active_to_slow / d5prime / myDF2$activesoil
soiltext2 <- (0.85 -  (0.996 - pas)) / 0.68
finesoil2 <- 0.75 - (1 - soiltext2)
finesoil2

# find psa coefficient; original value = 0.42
psa <- myDF2$slow_to_active / d6prime / myDF2$slowsoil

# find pap coefficient; original value = 0.004
pap <- myDF2$active_to_passive / d5prime / myDF2$activesoil

# find psp coefficient; original value = 0.03
psp <- myDF2$slow_to_passive / d6prime / myDF2$slowsoil

# find ppa coefficient; original value = 0.45
ppa <- myDF2$passive_to_active / d7prime / myDF2$passivesoil

# find kdec6 coefficient
factive <- myDF2$active_to_slow + myDF2$active_to_passive + myDF2$co2_rel_from_active_pool + myDF2$co2_released_exud
prime_y <- 0.6
prime_z <- 0.5
rt_slow_pool = (1.0 / prime_y) / (factive / (factive + prime_z))  # pmax(0.3, (factive / (factive + prime_z)))
kdec6_new <- 1 / rt_slow_pool
summary(kdec6_new)

# assign priming on coefficients
outDF[1,"prime_on"] <- pas[1]
outDF[2,"prime_on"] <- psa[1]
outDF[3,"prime_on"] <- pap[1]
outDF[4,"prime_on"] <- ppa[1]
outDF[5,"prime_on"] <- psp[1]
outDF[6,"prime_on"] <- mean(kdec6_new)


outDF


### Find relationship between extra C and kdec7
myDF2 <- fread("GDAY/outputs/Run9/Quasi_equil_model_spinup_equilib.csv",skip=1)
names(myDF2)

# compute extra C at each time step
myDF2$ar <- myDF2$cproot / (myDF2$cpleaf + myDF2$cproot + myDF2$cpstem)
myDF2$af <- myDF2$cpleaf / (myDF2$cpleaf + myDF2$cproot + myDF2$cpstem)
myDF2$aw <- myDF2$cpstem / (myDF2$cpleaf + myDF2$cproot + myDF2$cpstem)
myDF2$ariz <- 0.05 + 0.2 * pmax((myDF2$shoot/myDF2$shootn - 25)/25, 0)

myDF2$c_into_exud <- myDF2$npp * myDF2$ar * myDF2$ariz
myDF2$factive <- myDF2$active_to_slow + myDF2$active_to_passive + myDF2$co2_rel_from_active_pool + myDF2$co2_released_exud
myDF2$rt_passive <- 1/0.6/pmax(0.01, myDF2$factive/(myDF2$factive + 0.5))
myDF2$kdec7_pred <- (1 / myDF2$rt_passive) * 365

with(myDF2, plot(kdec7~c_into_exud))

with(myDF2, plot(kdec7_pred~kdec7))
