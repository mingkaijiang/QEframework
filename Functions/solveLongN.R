# Find the long term equilibrium nf and NPP under standard conditions - by finding the root
solveLongN <- function(CO2,Cpass,NinL) {
    fn <- function(nf) {
        solveNC(nf,allocn(nf)$af,CO2) - NConsLong(nf,allocn(nf),Cpass=Cpass,NinL)$NPP
    }
    equilnf <- uniroot(fn,interval=c(0.01,0.05))$root
    equilNPP <- solveNC(equilnf,af=allocn(equilnf)$af, CO2)
    ans <- data.frame(equilnf,equilNPP)
    return(ans)
}


# Find the long term equilibrium nf and NPP under standard conditions - by finding the root
# Specifically for respiration related calculations
solveLongN_respiration <- function(co2=350,Cpass,Nin, nwvar=T, nw=0.0005) {
    fn <- function(nf) {
        solveNC_respiration(nf,allocn(nf,nwvar=nwvar, nwood = nw),co2=co2) - NConsLong(nf,allocn(nf,nwvar=nwvar, nwood = nw),Cpass=Cpass,Nin=Nin)$NPP
    }
    equilnf <- uniroot(fn,interval=c(0.01,0.05))$root
    equilNPP <- solveNC_respiration(equilnf,adf=allocn(equilnf,nwvar=nwvar, nwood = nw), co2=co2)
    ans <- data.frame(equilnf,equilNPP)
    return(ans)
}

# Find the long term equilibrium nf and NPP under standard conditions - by finding the root
# specifically for CWD related calculuations
solveLongN_CWD <- function(co2=350,Cpass,Nin, nwvar=T) {
    fn <- function(nf) {
        solveNC(nf,allocn(nf,nwvar=nwvar)$af,co2=co2) - NConsLong_CWD(nf,allocn(nf,nwvar=nwvar),Cpass=Cpass,Nin=Nin)$NPP
    }
    equilnf <- uniroot(fn,interval=c(0.01,0.05))$root
    equilNPP <- solveNC(equilnf,af=allocn(equilnf,nwvar=nwvar)$af, co2=co2)
    ans <- data.frame(equilnf,equilNPP)
    return(ans)
}

# Find the long term equilibrium nf and NPP under standard conditions - by finding the root
# specifically for explicit mineral pools
solveLongN_expl_min <- function(co2=350,Cpass,Nin, nwvar=T) {
    fn <- function(nf) {
        solveNC(nf,allocn(nf,nwvar=nwvar)$af,co2=co2) - NConsLong_expl_min(nf,allocn(nf,nwvar=nwvar),Cpass=Cpass,Nin=Nin)$NPP
    }
    equilnf <- uniroot(fn,interval=c(0.01,0.05))$root
    equilNPP <- solveNC(equilnf,af=allocn(equilnf,nwvar=nwvar)$af, co2=co2)
    ans <- data.frame(equilnf,equilNPP)
    return(ans)
}

# Find the long term equilibrium nf and NPP under standard conditions - by finding the root
# specifically for nuptake ~ root biomass  - O-CN approach
# i.e. N uptake as a saturating function of mineral N
solveLongN_root_ocn <- function(co2=350, Nin, nwvar=F) {
    fn <- function(nf) {
        solveNC(nf,allocn(nf,nwvar=nwvar)$af,co2=co2) - NConsLong_root_ocn(nf,allocn(nf,nwvar=nwvar),Nin=Nin)$NPP
    }
    equilnf <- uniroot(fn,interval=c(0.01,0.05))$root
    equilNPP <- solveNC(equilnf,af=allocn(equilnf,nwvar=nwvar)$af, co2=co2)
    ans <- data.frame(equilnf,equilNPP)
    return(ans)
}


# Find the long term equilibrium nf and NPP under standard conditions - by finding the root
# specifically for nuptake ~ root biomass  - O-CN approach
# i.e. N uptake as a saturating function of mineral N
solveLongN_root_gday <- function(co2=350,Cpass,Nin, nwvar=T) {
    fn <- function(nf) {
        solveNC(nf,allocn(nf,nwvar=nwvar)$af,co2=co2) - NConsLong_root_gday(nf,allocn(nf,nwvar=nwvar),Cpass=Cpass,Nin=Nin)$NPP
    }
    equilnf <- uniroot(fn,interval=c(0.01,0.05))$root
    equilNPP <- solveNC(equilnf,af=allocn(equilnf,nwvar=nwvar)$af, co2=co2)
    ans <- data.frame(equilnf,equilNPP)
    return(ans)
}

# Find the long term equilibrium nf and NPP under standard conditions - by finding the root
# specifically for variable passive pool stoichiometry
solveLongN_variable_pass <- function(co2=350,Cpass,Nin, nwvar=T, eqNPP) {
    fn <- function(nf) {
        solveNC(nf,allocn(nf,nwvar=nwvar)$af,co2=co2) - NConsLong_variable_pass(nf,allocn(nf,nwvar=nwvar),Cpass=Cpass,Nin=Nin, eqNPP=eqNPP)$NPP
    }
    equilnf <- uniroot(fn,interval=c(0.01,0.05))$root
    equilNPP <- solveNC(equilnf,af=allocn(equilnf,nwvar=nwvar)$af, co2=co2)
    ans <- data.frame(equilnf,equilNPP)
    return(ans)
}


# Find the long term equilibrium nf and NPP under standard conditions - by finding the root
# specifically for exudation
solveLongN_exudation <- function(co2=350,Cpass,Nin, nwvar=T) {
    fn <- function(nf) {
        solveNC(nf,allocn_exudation(nf,nwvar=nwvar)$af,co2=co2) - NConsLong_exudation(nf,allocn_exudation(nf,nwvar=nwvar),Cpass=Cpass,Nin=Nin)$NPP
    }
    equilnf <- uniroot(fn,interval=c(0.01,0.05))$root
    equilNPP <- solveNC(equilnf,af=allocn_exudation(equilnf,nwvar=nwvar)$af, co2=co2)
    ans <- data.frame(equilnf,equilNPP)
    return(ans)
}