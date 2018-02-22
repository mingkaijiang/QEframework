VL_constraint_baseline_CLM <- function(df,a) {
    ### parameters:
    ###            f: fraction of plant enters into soil
    ###            nsoil: soil NC ratio

    # equation for N constraint with just leaching
    U0 <- Nin
    
    # Calculate Nmin
    Nmin <- Nin / leachn
    
    # Calculate NPP
    NPP_act <- Nmin / ((1-f) * (a$nfl*a$af + a$nr*a$ar + a$nw*a$aw) + nsoil)
 
    # return in kg C m-2 yr-1
    NPP <- NPP_act *10^-3    
    
    browser()
    
    # out df
    out <- data.frame(NPP,Nmin)
    colnames(out) <- c("NPP", "nmin")
    
    return(out)   
}