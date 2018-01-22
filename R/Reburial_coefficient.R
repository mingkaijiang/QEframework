### Function to look at priming effect on reburial coefficient

reburial_comparison <- function() {
    
    #### Baseline case - basic N model
    source("Parameters/Analytical_Run2_Parameters.R")
    nfseq <- seq(0.001,0.1,by=0.001)
    a_vec <- allocn(nfseq)
    PC350 <- photo_constraint_full_cn(nfseq,a_vec,CO2=CO2_1)
    PC700 <- photo_constraint_full_cn(nfseq,a_vec,CO2=CO2_2)
    NCVLONG <- VLong_constraint_N(nfseq,a_vec)
    VLong <- solveVLong_full_cn(CO2=CO2_1)
    aequil <- allocn(VLong$equilnf)
    pass_case1 <- slow_pool(df=VLong$equilnf, a=aequil)
    omegap_case1 <- aequil$af*pass_case1$omegafp + aequil$ar*pass_case1$omegarp
    omegas_case1 <- aequil$af*pass_case1$omegafs + aequil$ar*pass_case1$omegars
    CpassVLong_case1 <- omegap_case1*VLong$equilNPP/pass_case1$decomp_p/(1-pass_case1$qpq)*1000.0
    NCHUGH <- Long_constraint_N(nfseq,a_vec,Cpass=CpassVLong_case1,Nin)
    equil_long_350 <- solveLong_full_cn(CO2=CO2_1, Cpass=CpassVLong_case1, NinL = Nin)
    CslowLong_case1 <- omegas_case1*equil_long_350$equilNPP/pass_case1$decomp_s/(1-pass_case1$qsq)*1000.0
    
    #### Priming turned off, exudation turned on
    source("Parameters/Analytical_Run10_Parameters.R")
    nfseq <- round(seq(0.001, 0.1, by = 0.001),5)
    a_vec <- as.data.frame(allocn_exudation(nfseq))
    PC350 <- photo_constraint_full_cn(nfseq,a_vec,CO2=CO2_1)
    PC700 <- photo_constraint_full_cn(nfseq,a_vec,CO2=CO2_2)
    NCVLONG <- NConsVLong(df=nfseq,a=a_vec)
    VLongN <- solveVLong_full_cn_medium(CO2_1)
    aequiln <- allocn_exudation(VLongN$equilnf)
    pass_case2 <- slow_pool(df=VLongN$equilnf, a=aequiln)
    omegap_case2 <- aequiln$af*pass_case2$omegafp + aequiln$ar*pass_case2$omegarp
    omegas_case2 <- aequiln$af*pass_case2$omegafs + aequiln$ar*pass_case2$omegars
    CpassVLong_case2 <- omegap_case2*VLongN$equilNPP/pass_case2$decomp_p/(1-pass_case2$qpq)*1000.0
    CslowLong_case2 <- omegas_case2*VLongN$equilNPP/pass_case2$decomp_s/(1-pass_case2$qsq)*1000.0
    
    
    #### Priming turned on, exudation turned on
    source("Parameters/Analytical_Run10_Parameters.R")
    nfseq <- round(seq(0.001, 0.1, by = 0.001),5)
    a_vec <- as.data.frame(allocn_exudation(nfseq))
    PC350 <- photo_constraint_full_cn(nfseq,a_vec,CO2=CO2_1)
    PC700 <- photo_constraint_full_cn(nfseq,a_vec,CO2=CO2_2)
    NCVLONG <- NConsVLong(df=nfseq,a=a_vec)
    VLongN <- solveVLong_full_cn_medium(CO2_1)
    aequiln <- allocn_exudation(VLongN$equilnf)
    pass_case3 <- slow_pool(df=VLongN$equilnf, a=aequiln)
    omegap_case3 <- aequiln$af*pass_case3$omegafp + aequiln$ar*pass_case3$omegarp
    omegas_case3 <- aequiln$af*pass_case3$omegafs + aequiln$ar*pass_case3$omegars
    CpassVLong_case3 <- omegap_case3*VLongN$equilNPP/pass_case3$decomp_p/(1-pass_case3$qpq)*1000.0
    CslowLong_no_priming <- omegas_case3*VLongN$equilNPP/pass_case3$decomp_s/(1-pass_case3$qsq)*1000.0
    NrelwoodVLong <- aequiln$aw*aequiln$nw*VLongN$equilNPP*1000.0
    NCHUGH <- NConsLong(nfseq, a_vec, CpassVLong_case3,
                        NinL = Nin)
    c_into_active <- VLongN$equilNPP * aequiln$ar * aequiln$ariz * rhizo_cue * 1000.0
    n_into_active <- c_into_active * aequiln$nr
    n_active_gap <- c_into_active * nca - n_into_active
    new_kdec <- pass_case3$decomp_s * (1 + km) * pmax(c_into_active/(c_into_active + km), 0.3)
    CslowLong_case3 <- omegas_case3*VLongN$equilNPP/new_kdec/(1-pass_case3$qsq)*1000.0
    
    #### save output into a table
    out.table <- rbind(pass_case1, pass_case2, pass_case3)
    out.table <- as.data.frame(out.table)
    out.table$case <- c("baseline", "exudation", "priming")
    out.table$decomp_s[3] <- new_kdec
    out.table$omegap <- c(omegap_case1, omegap_case2, omegap_case3)
    out.table$omegas <- c(omegas_case1, omegas_case2, omegas_case3)
    out.table$Cpass <- c(CpassVLong_case1, CpassVLong_case2, CpassVLong_case3)
    out.table$Cslow <- c(CslowLong_case1, CslowLong_case2, CslowLong_case3)
    
    #### organize the out table
    out.table <- out.table[,c("case", "decomp_p", "decomp_s", "qpq", "qsq",
                             "omegafp", "omegarp", "omegafs", "omegars", 
                             "omegap", "omegas", "Cpass", "Cslow")]
    
    #### drop case2
    out.table <- out.table[-2,]
    
    #### save table as supplementary table
    write.csv(out.table, "Tables/TableS3.csv", row.names=F)
}

### Script
reburial_comparison()
