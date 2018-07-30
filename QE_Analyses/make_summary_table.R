######################## Functions ###################################
make_summary_table <- function(destDir) {
    
    if(!dir.exists(destDir)) {
        dir.create(destDir)
    }
    
    ### Make the output table
    out.tab <- matrix(ncol=7, nrow = 12)
    colnames(out.tab) <- c("model", "NPP_350", "NPP_700", "I", "M", "L", "VL")
    out.tab <- as.data.frame(out.tab)
    out.tab[1,"model"] <- "Baseline variable wood"
    out.tab[2,"model"] <- "Baseline fixed wood"
    out.tab[3,"model"] <- "Explicit N uptake, fixed coefficient, variable wood"
    out.tab[4,"model"] <- "Explicit N uptake, fixed coefficient, fixed wood"
    out.tab[5,"model"] <- "Explicit N uptake, GDAY, variable wood"
    out.tab[6,"model"] <- "Explicit N uptake, OCN, variable wood"
    out.tab[7,"model"] <- "Priming, variable wood"
    out.tab[8,"model"] <- "Priming, fixed wood"
    out.tab[9,"model"] <- "Potential NPP, variable wood"
    out.tab[10,"model"] <- "Potential NPP, fixed wood"
    #out.tab[11,"model"] <- "C cost of N uptake, variable wood"
    #out.tab[12,"model"] <- "C cost of N uptake, fixed wood"
    
    
    ### Run 1 Baseline variable wood
    out.tab[1,"NPP_350"] <- round(r1$eDF$NPP_VL[1],2)
    out.tab[1,"NPP_700"] <-  round(r1$eDF$NPP_VL[2],2)
    out.tab[1,"I"] <- round((r1$eDF$NPP_I[2] - r1$eDF$NPP_VL[1]) / r1$eDF$NPP_VL[1] * 100,2)
    out.tab[1,"M"] <- round((r1$eDF$NPP_M[2] - r1$eDF$NPP_VL[1]) / r1$eDF$NPP_VL[1] * 100,2)
    out.tab[1,"L"] <- round((r1$eDF$NPP_L[2] - r1$eDF$NPP_VL[1]) / r1$eDF$NPP_VL[1] * 100,2)
    out.tab[1,"VL"] <- round((r1$eDF$NPP_VL[2] - r1$eDF$NPP_VL[1]) / r1$eDF$NPP_VL[1] * 100,2)
    
    ### Run 2 Baseline fixed wood
    out.tab[2,"NPP_350"] <- round(r2$eDF$NPP_VL[1],2)
    out.tab[2,"NPP_700"] <-  round(r2$eDF$NPP_VL[2],2)
    out.tab[2,"I"] <- round((r2$eDF$NPP_I[2] - r2$eDF$NPP_VL[1]) / r2$eDF$NPP_VL[1] * 100,2)
    out.tab[2,"M"] <- round((r2$eDF$NPP_M[2] - r2$eDF$NPP_VL[1]) / r2$eDF$NPP_VL[1] * 100,2)
    out.tab[2,"L"] <- round((r2$eDF$NPP_L[2] - r2$eDF$NPP_VL[1]) / r2$eDF$NPP_VL[1] * 100,2)
    out.tab[2,"VL"] <- round((r2$eDF$NPP_VL[2] - r2$eDF$NPP_VL[1]) / r2$eDF$NPP_VL[1] * 100,2)
    
    ### Run 3 Explicit N uptake, fixed coefficient, variable wood
    out.tab[3,"NPP_350"] <- round(r3$eDF$NPP_VL[1],2)
    out.tab[3,"NPP_700"] <-  round(r3$eDF$NPP_VL[2],2)
    out.tab[3,"I"] <- round((r3$eDF$NPP_I[2] - r3$eDF$NPP_VL[1]) / r3$eDF$NPP_VL[1] * 100,2)
    out.tab[3,"M"] <- round((r3$eDF$NPP_M[2] - r3$eDF$NPP_VL[1]) / r3$eDF$NPP_VL[1] * 100,2)
    out.tab[3,"L"] <- round((r3$eDF$NPP_L[2] - r3$eDF$NPP_VL[1]) / r3$eDF$NPP_VL[1] * 100,2)
    out.tab[3,"VL"] <- round((r3$eDF$NPP_VL[2] - r3$eDF$NPP_VL[1]) / r3$eDF$NPP_VL[1] * 100,2)
    
    ### Run 4 Explicit N uptake, fixed coefficient, fixed wood
    out.tab[4,"NPP_350"] <- round(r4$eDF$NPP_VL[1],2)
    out.tab[4,"NPP_700"] <-  round(r4$eDF$NPP_VL[2],2)
    out.tab[4,"I"] <- round((r4$eDF$NPP_I[2] - r4$eDF$NPP_VL[1]) / r4$eDF$NPP_VL[1] * 100,2)
    out.tab[4,"M"] <- round((r4$eDF$NPP_M[2] - r4$eDF$NPP_VL[1]) / r4$eDF$NPP_VL[1] * 100,2)
    out.tab[4,"L"] <- round((r4$eDF$NPP_L[2] - r4$eDF$NPP_VL[1]) / r4$eDF$NPP_VL[1] * 100,2)
    out.tab[4,"VL"] <- round((r4$eDF$NPP_VL[2] - r4$eDF$NPP_VL[1]) / r4$eDF$NPP_VL[1] * 100,2)
    
    ### Run 5 Explicit N uptake, GDAY, variable wood
    out.tab[5,"NPP_350"] <- round(r5$eDF$NPP_VL[1],2)
    out.tab[5,"NPP_700"] <-  round(r5$eDF$NPP_VL[2],2)
    out.tab[5,"I"] <- round((r5$eDF$NPP_I[2] - r5$eDF$NPP_VL[1]) / r5$eDF$NPP_VL[1] * 100,2)
    out.tab[5,"M"] <- round((r5$eDF$NPP_M[2] - r5$eDF$NPP_VL[1]) / r5$eDF$NPP_VL[1] * 100,2)
    out.tab[5,"L"] <- round((r5$eDF$NPP_L[2] - r5$eDF$NPP_VL[1]) / r5$eDF$NPP_VL[1] * 100,2)
    out.tab[5,"VL"] <- round((r5$eDF$NPP_VL[2] - r5$eDF$NPP_VL[1]) / r5$eDF$NPP_VL[1] * 100,2)
    
    ### Run 6 Explicit N uptake, OCN, variable wood
    out.tab[6,"NPP_350"] <- round(r6$eDF$NPP_VL[1],2)
    out.tab[6,"NPP_700"] <-  round(r6$eDF$NPP_VL[2],2)
    out.tab[6,"I"] <- round((r6$eDF$NPP_I[2] - r6$eDF$NPP_VL[1]) / r6$eDF$NPP_VL[1] * 100,2)
    out.tab[6,"M"] <- round((r6$eDF$NPP_M[2] - r6$eDF$NPP_VL[1]) / r6$eDF$NPP_VL[1] * 100,2)
    out.tab[6,"L"] <- round((r6$eDF$NPP_L[2] - r6$eDF$NPP_VL[1]) / r6$eDF$NPP_VL[1] * 100,2)
    out.tab[6,"VL"] <- round((r6$eDF$NPP_VL[2] - r6$eDF$NPP_VL[1]) / r6$eDF$NPP_VL[1] * 100,2)
    
    ### Run 7 Priming, fixed wood
    out.tab[8,"NPP_350"] <- round(r7$eDF$NPP_VL[1],2)
    out.tab[8,"NPP_700"] <-  round(r7$eDF$NPP_VL[2],2)
    out.tab[8,"I"] <- round((r7$eDF$NPP_I[2] - r7$eDF$NPP_VL[1]) / r7$eDF$NPP_VL[1] * 100,2)
    out.tab[8,"M"] <- round((r7$eDF$NPP_M[2] - r7$eDF$NPP_VL[1]) / r7$eDF$NPP_VL[1] * 100,2)
    out.tab[8,"L"] <- round((r7$eDF$NPP_L[2] - r7$eDF$NPP_VL[1]) / r7$eDF$NPP_VL[1] * 100,2)
    out.tab[8,"VL"] <- round((r7$eDF$NPP_VL[2] - r7$eDF$NPP_VL[1]) / r7$eDF$NPP_VL[1] * 100,2)
    
    ### Run 8 Priming, variable wood
    out.tab[7,"NPP_350"] <- round(r8$eDF$NPP_VL[1],2)
    out.tab[7,"NPP_700"] <-  round(r8$eDF$NPP_VL[2],2)
    out.tab[7,"I"] <- round((r8$eDF$NPP_I[2] - r8$eDF$NPP_VL[1]) / r8$eDF$NPP_VL[1] * 100,2)
    out.tab[7,"M"] <- round((r8$eDF$NPP_M[2] - r8$eDF$NPP_VL[1]) / r8$eDF$NPP_VL[1] * 100,2)
    out.tab[7,"L"] <- round((r8$eDF$NPP_L[2] - r8$eDF$NPP_VL[1]) / r8$eDF$NPP_VL[1] * 100,2)
    out.tab[7,"VL"] <- round((r8$eDF$NPP_VL[2] - r8$eDF$NPP_VL[1]) / r8$eDF$NPP_VL[1] * 100,2)
    
    ### Run 9 Potential NPP, fixed wood
    out.tab[10,"NPP_350"] <- round(r9$eDF$NPP_act_VL[1],2)
    out.tab[10,"NPP_700"] <-  round(r9$eDF$NPP_act_VL[2],2)
    out.tab[10,"I"] <- round((r9$eDF$NPP_I[2] - r9$eDF$NPP_act_VL[1]) / r9$eDF$NPP_act_VL[1] * 100,2)
    out.tab[10,"M"] <- round((r9$eDF$NPP_act_M[2] - r9$eDF$NPP_act_VL[1]) / r9$eDF$NPP_act_VL[1] * 100,2)
    out.tab[10,"L"] <- round((r9$eDF$NPP_act_L[2] - r9$eDF$NPP_act_VL[1]) / r9$eDF$NPP_act_VL[1] * 100,2)
    out.tab[10,"VL"] <- round((r9$eDF$NPP_act_VL[2] - r9$eDF$NPP_act_VL[1]) / r9$eDF$NPP_act_VL[1] * 100,2)
    
    ### Run 10 Potential NPP, variable wood
    out.tab[9,"NPP_350"] <- round(r10$eDF$NPP_act_VL[1],2)
    out.tab[9,"NPP_700"] <-  round(r10$eDF$NPP_act_VL[2],2)
    out.tab[9,"I"] <- round((r10$eDF$NPP_I[2] - r10$eDF$NPP_act_VL[1]) / r10$eDF$NPP_act_VL[1] * 100,2)
    out.tab[9,"M"] <- round((r10$eDF$NPP_act_M[2] - r10$eDF$NPP_act_VL[1]) / r10$eDF$NPP_act_VL[1] * 100,2)
    out.tab[9,"L"] <- round((r10$eDF$NPP_act_L[2] - r10$eDF$NPP_act_VL[1]) / r10$eDF$NPP_act_VL[1] * 100,2)
    out.tab[9,"VL"] <- round((r10$eDF$NPP_act_VL[2] - r10$eDF$NPP_act_VL[1]) / r10$eDF$NPP_act_VL[1] * 100,2)
    
    #### Run 11 C cost of N uptake, fixed wood
    #out.tab[12,"NPP_350"] <- round(r11$eDF$NPP_VL[1],2)
    #out.tab[12,"NPP_700"] <-  round(r11$eDF$NPP_VL[2],2)
    #out.tab[12,"I"] <- round((r11$eDF$NPP_I[2] - r11$eDF$NPP_VL[1]) / r11$eDF$NPP_VL[1] * 100,2)
    #out.tab[12,"M"] <- round((r11$eDF$NPP_M[2] - r11$eDF$NPP_VL[1]) / r11$eDF$NPP_VL[1] * 100,2)
    #out.tab[12,"L"] <- round((r11$eDF$NPP_L[2] - r11$eDF$NPP_VL[1]) / r11$eDF$NPP_VL[1] * 100,2)
    #out.tab[12,"VL"] <- round((r11$eDF$NPP_VL[2] - r11$eDF$NPP_VL[1]) / r11$eDF$NPP_VL[1] * 100,2)
    #
    #### Run 12 C cost of N uptake, variable
    #out.tab[11,"NPP_350"] <- round(r12$eDF$NPP_VL[1],2)
    #out.tab[11,"NPP_700"] <-  round(r12$eDF$NPP_VL[2],2)
    #out.tab[11,"I"] <- round((r12$eDF$NPP_I[2] - r12$eDF$NPP_VL[1]) / r12$eDF$NPP_VL[1] * 100,2)
    #out.tab[11,"M"] <- round((r12$eDF$NPP_M[2] - r12$eDF$NPP_VL[1]) / r12$eDF$NPP_VL[1] * 100,2)
    #out.tab[11,"L"] <- round((r12$eDF$NPP_L[2] - r12$eDF$NPP_VL[1]) / r12$eDF$NPP_VL[1] * 100,2)
    #out.tab[11,"VL"] <- round((r12$eDF$NPP_VL[2] - r12$eDF$NPP_VL[1]) / r12$eDF$NPP_VL[1] * 100,2)
    
    ### Save the output table
    write.csv(out.tab, paste0(destDir, "/summary_table.csv"), row.names=F)
    
    return(out.tab)
}