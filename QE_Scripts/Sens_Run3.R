
#### Sensitivity of Run 3 leaching 
#### Plot NPP vs. Nleaching at different timescales
####
#### Assumptions:
#### 1. Variable wood NC
#### 2. baseline N cycle
#### 3. Explicit mineral N pool
################################################################################
#### Functions
Sens_Run3 <- function() {
    
    ### Get 
    r3 <- Perform_Analytical_Run3(f.flag = 2)
    
    source("Parameters/Analytical_Run3_Parameters.R")
    
    
    ### split the data.list into two data frames
    cdf <- r3$cDF
    edf <- r3$eDF
    
    ### extract CO2 conc. from constraint DF
    cdf1 <- cdf[cdf$CO2 == 400, ]
    cdf2 <- cdf[cdf$CO2 == 800, ]
    
    ### extract CO2 conc. from eq DF
    edf1 <- edf[edf$CO2 == 400, ]
    edf2 <- edf[edf$CO2 == 800, ]
    
    ### Create a plot df
    plotDF <- data.frame(c("VL_400", "M", "L", "VL_800"), NA, NA)
    colnames(plotDF) <- c("Timeframe", "NPP", "Nleach")
    
    ### Add VL leaching rate
    plotDF$NPP[plotDF$Timeframe=="VL_400"] <- edf1$NPP_VL
    plotDF$NPP[plotDF$Timeframe=="VL_800"] <- edf2$NPP_VL
    plotDF$Nleach[plotDF$Timeframe=="VL_400"] <- Nin 
    plotDF$Nleach[plotDF$Timeframe=="VL_800"] <- Nin
    plotDF$NPP[plotDF$Timeframe=="L"] <- edf2$NPP_L
    plotDF$NPP[plotDF$Timeframe=="M"] <- edf2$NPP_M
    
    ### work out nleaching L
    nfseq <- round(seq(0.001, 0.1, by = 0.001),5)
    a_nf <- as.data.frame(alloc(nfseq))
    P350 <- photo_constraint_full(nf=nfseq, nfdf=a_nf, CO2=CO2_1)
    VL <- VL_constraint_expl_min(df=nfseq, a=a_nf)
    VL_eq <- solve_VL_full_expl_min(CO2=CO2_1)
    a_eq <- alloc(VL_eq$nf)
    s_coef <- soil_coef(df=VL_eq$nf, a=a_eq)
    omega_ap <- a_eq$af*s_coef$omega_af_pass + a_eq$ar*s_coef$omega_ar_pass + a_eq$aw*s_coef$omega_aw_pass
    omega_as <- a_eq$af*s_coef$omega_af_slow + a_eq$ar*s_coef$omega_ar_slow + a_eq$aw*s_coef$omega_aw_slow
    C_pass_VL <- omega_ap*VL_eq$NPP/s_coef$decomp_pass/(1-s_coef$qq_pass)*1000.0
    L <- L_constraint_expl_min(df=nfseq, a=a_nf, 
                               C_pass=C_pass_VL, Nin_L = Nin)

    test <- round(edf2$NPP_L, 1)
    loc1 <- 17
    loc2 <- 18
    mydf <- L[loc1:loc2, c("NPP", "nleach")]
    lm.model <- lm(mydf[,2]~mydf[,1])

    # calculate 
    a <- as.numeric(coefficients(lm.model)[1])
    b <- as.numeric(coefficients(lm.model)[2])
    
    eq_Nleach_L <- a+b*edf2$NPP_L
    plotDF$Nleach[plotDF$Timeframe=="L"] <- eq_Nleach_L
    
    ### obtain M line
    L_eq <- solve_L_full_expl_min(CO2=CO2_1, C_pass=C_pass_VL, 
                                  Nin_L = Nin)
    C_slow_L <- omega_as*L_eq$NPP/s_coef$decomp_slow/(1-s_coef$qq_slow)*1000.0
    N_wood_L <- a_eq$aw*a_eq$nw*VL_eq$NPP*1000.0
    M <- M_constraint_expl_min(df=nfseq,a=a_nf, 
                               C_pass=C_pass_VL, 
                               C_slow=C_slow_L, 
                               Nin_L = Nin+N_wood_L)
        
    test <- round(edf2$NPP_M, 1)
    loc1 <- 9
    loc2 <- 10
    mydf <- M[loc1:loc2, c("NPP", "nleach")]
    lm.model <- lm(mydf[,2]~mydf[,1])
    
    # calculate 
    a <- as.numeric(coefficients(lm.model)[1])
    b <- as.numeric(coefficients(lm.model)[2])
    
    eq_Nleach_M <- a+b*edf2$NPP_M
    plotDF$Nleach[plotDF$Timeframe=="M"] <- eq_Nleach_M
    
    ### Plotting
    cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
    
    p<-ggplot() + 
        geom_bar(data=plotDF, aes(x=as.factor(NPP), y=Nleach, fill=Timeframe), 
                 stat="identity", position=position_dodge()) +
        theme_linedraw() +
        theme(panel.grid.minor=element_blank(),
              axis.text=element_text(size=14),
              axis.title=element_text(size=16),
              legend.text=element_text(size=14),
              legend.title=element_text(size=16),
              axis.text.x = element_text(angle = 90, hjust=1),
              panel.grid.major=element_line(color="grey")) + 
        ylab(expression(paste("N leaching rate (kg N ", m^-2, " ", yr^-1, ")")))+ xlab("") + 
        scale_fill_manual(name="Timescale", 
                          values = c("VL_400" = cbPalette[2], "M" = cbPalette[3],
                                      "L" = cbPalette[4], "VL_800" = cbPalette[7]),
                          labels = c(expression(VL[400]), "M", "L", expression(VL[800]))) 
    plot(p)

    pdf("Plots/sensitivity_NPP_Nleach.pdf", width=10, height=6)
    plot(p)
    dev.off()
    
    write.csv(plotDF, "Tables/sens_NPP_Nleach.csv")
    
}

