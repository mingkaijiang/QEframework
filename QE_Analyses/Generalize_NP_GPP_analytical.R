######### Obtainning relationship between leaf NC and leaf PC and LUE for N and NP models

source("R/prepare_R.R")

#### NP model
source("Parameters/Analytical_Run1_Parameters.R")

# create a range of nc for shoot to initiate
nfseq <- round(seq(0.01, 0.1, by = 0.001),5)
a_nf <- as.data.frame(allocn(nfseq))

# using very long term relationship to calculate pf from nf
pfseq <- inferpfVL(nfseq, a_nf)
a_pf <- as.data.frame(allocp(pfseq))

# create storage df
outDF <- data.frame(rep(nfseq,10), rep(pfseq, 10), NA, NA)
colnames(outDF) <- c("nfseq", "pfseq", "LUE", "CO2")
co2.list <- seq(350, 800, by = 50)
outDF$CO2 <- rep(co2.list, each = 91)

# calculate photosynthetic constraint at CO2 = 350

for (i in 1:length(co2.list)) {
    Photo350 <- photo_constraint_full_cnp(nfseq, pfseq, a_nf, a_pf, co2.list[i])
    LUE_NP <- LUE_full_cnp_walker(nfseq, a_pf, pfseq, CO2=co2.list[i], Photo350*1000.0) 
    outDF[outDF$CO2 == co2.list[i], "LUE"] <- LUE_NP
}

with(outDF, plot(nfseq, LUE))
with(outDF, plot(pfseq, LUE))

ratio1 <- (outDF[outDF$CO2 == 400, "LUE"] - outDF[outDF$CO2 == 350, "LUE"]) / outDF[outDF$CO2 == 350, "LUE"]
ratio2 <- (outDF[outDF$CO2 == 450, "LUE"] - outDF[outDF$CO2 == 400, "LUE"]) / outDF[outDF$CO2 == 400, "LUE"]

# trial one
lm.np <- lm(log(outDF$LUE)~log(outDF$CO2) + log(outDF$nfseq) + (log(outDF$nfseq):log(outDF$pfseq)))
summary(lm.np)

#lue_np_pred <- 0.0288 + 0.82 * nfseq - 118.1 * nfseq * pfseq
i <- 10
lue_np_pred <- exp(-5.06 + 0.21 * log(co2.list[i]) - 1.24 * log(nfseq) - 0.15 * log(nfseq) * log(pfseq)) 
Photo350 <- photo_constraint_full_cnp(nfseq, pfseq, a_nf, a_pf, co2.list[i])
LUE_NP <- LUE_full_cnp_walker(nfseq, a_pf, pfseq, CO2=co2.list[i], Photo350*1000.0) 
plot(lue_np_pred~LUE_NP)
abline(a=0,b=1)

# trial two
lm.np <- lm(LUE_NP~-1+nfseq+nfseq:pfseq)
summary(lm.np)

lue_np_pred <- 1.92 * nfseq - 309.63 * nfseq * pfseq
plot(lue_np_pred~LUE_NP)
abline(a=0,b=1)

# trial three
lm.np <- lm(LUE_NP~nfseq:pfseq)
summary(lm.np)

lue_np_pred <- 0.047 +  34.02 * nfseq * pfseq
plot(lue_np_pred~LUE_NP)
abline(a=0,b=1)

# trial four
lm.np <- lm(LUE_NP~nfseq+pfseq+nfseq:pfseq)
summary(lm.np)

lue_np_pred <- 0.029 + 1946 * nfseq -41540 * pfseq - 117.8 * nfseq * pfseq
plot(lue_np_pred~LUE_NP)
abline(a=0,b=1)

# trial 5
lm.np <- lm(log(LUE_NP)~log(nfseq) + log(nfseq):log(pfseq))
summary(lm.np)

lue_np_pred <- exp(-3.85 - 1.25 * log(nfseq) - 0.15 * log(pfseq) * log(nfseq))
plot(lue_np_pred~LUE_NP)
abline(a=0,b=1)

# trial 6
lm.np <- lm(log(LUE_NP)~log(pfseq) + log(nfseq):log(pfseq))
summary(lm.np)

lue_np_pred <- exp(-7.69 - 1.25 * log(pfseq) - 0.15 * log(pfseq) * log(nfseq))
plot(lue_np_pred~LUE_NP)
abline(a=0,b=1)

# trial 7
lm.np <- lm(log(LUE_NP)~log(nfseq):log(pfseq))
summary(lm.np)

lue_np_pred <- exp(-2.49 - 0.02 * log(pfseq) * log(nfseq))
plot(lue_np_pred~LUE_NP)
abline(a=0,b=1)


# trial 8
lm.np <- lm(log(LUE_NP)~log(pfseq) + log(nfseq))
summary(lm.np)

lue_np_pred <- exp(-17588 - 5745 * log(pfseq) + 5745 * log(nfseq))
plot(lue_np_pred~LUE_NP)
abline(a=0,b=1)

# best model: trial 1 or trial 5

#### NP model
source("Parameters/Analytical_Run2_Parameters.R")

# create a range of nc for shoot to initiate
nfseq <- round(seq(0.01, 0.1, by = 0.001),5)
a_nf <- as.data.frame(allocn(nfseq))

# create storage df
outDF <- data.frame(rep(nfseq,10), NA, NA)
colnames(outDF) <- c("nfseq", "LUE", "CO2")
co2.list <- seq(350, 800, by = 50)
outDF$CO2 <- rep(co2.list, each = 91)

# calculate photosynthetic constraint at CO2 = 350

for (i in 1:length(co2.list)) {
    Photo350 <- photo_constraint_full_cn(nfseq,  a_nf, co2.list[i])
    LUE_N <- LUE_full_cn_walker(nfseq, a_nf, CO2=co2.list[i], Photo350*1000.0) 
    outDF[outDF$CO2 == co2.list[i], "LUE"] <- LUE_N
}

# trial one
lm.n <- lm(outDF$LUE~outDF$nfseq)
summary(lm.n)

lue_n_pred <- 0.04 + 0.09 * nfseq
plot(lue_n_pred~LUE_N)
abline(a=0,b=1)

# trial two
lm.n <- lm(log(outDF$LUE)~log(outDF$nfseq))
summary(lm.n)

lue_n_pred <- exp(-2.76 + 0.09 * log(nfseq))
plot(lue_n_pred~LUE_N)
abline(a=0,b=1)

# trial three
lm.n <- lm(log(outDF$LUE)~log(outDF$CO2) + log(outDF$nfseq) + (log(outDF$nfseq):log(outDF$CO2)))
summary(lm.n)

lue_n_pred <- exp(-4.03 + 0.22 * log(co2.list[i]) + 0.14 * log(nfseq) - 0.01 * log(nfseq)*log(co2.list[i]))
plot(lue_n_pred~LUE_N)
abline(a=0,b=1)

# trial four
lm.n <- lm(log(outDF$LUE)~log(outDF$CO2) + log(outDF$nfseq))
summary(lm.n)

lue_n_pred <- exp(-4.23 + 0.25 * log(co2.list[i]) + 0.08 * log(nfseq))
plot(lue_n_pred~LUE_N)
abline(a=0,b=1)

# trial five
lm.n <- lm(log(outDF$LUE)~log(outDF$CO2) + (log(outDF$nfseq):log(outDF$CO2)))
summary(lm.n)

lue_n_pred <- exp(-4.47 + 0.29 * log(co2.list[i]) + 0.01 * log(nfseq) * log(co2.list[i]))
plot(lue_n_pred~LUE_N)
abline(a=0,b=1)


# trial six
lm.n <- lm(log(outDF$LUE)~(log(outDF$nfseq):log(outDF$CO2)))
summary(lm.n)

lue_n_pred <- exp(-2.74 + 0.01 * log(nfseq) * log(co2.list[i]))
plot(lue_n_pred~LUE_N)
abline(a=0,b=1)


#### Save 1:1 line comparisons of the simple model vs. complex photosynthesis model
tiff("Plots/LUE_simple_vs_complex.tiff",
     width = 10, height = 5, units = "in", res = 300)
par(mfrow=c(1,2), mar=c(5.1,6.1,2.1,2.1))

# NP model
lm.np <- lm(log(LUE_NP)~log(nfseq) + log(nfseq):log(pfseq))
summary(lm.np)

lue_np_pred <- exp(-3.85 - 1.25 * log(nfseq) - 0.15 * log(pfseq) * log(nfseq))
plot(lue_np_pred~LUE_NP, xlab = "LUE obs", ylab = "LUE pred", type="p")
abline(a=0,b=1, col="red", lty=2)

# N model
lm.n <- lm(log(LUE_N)~log(nfseq))
summary(lm.n)

lue_n_pred <- exp(-2.76 + 0.09 * log(nfseq))
plot(lue_n_pred~LUE_N, xlab = "LUE obs", ylab = "LUE pred", type="p")
abline(a=0,b=1, col="red", lty=2)

dev.off()