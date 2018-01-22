calculate_michaelis_menten_parameter <- function(tk,mt) {
    
    # Michaelis-Menten coefficents for carboxylation by Rubisco 
    Kc = arrh(mt, kc25, eac, tk);
    
    # Michaelis-Menten coefficents for oxygenation by Rubisco 
    Ko = arrh(mt, ko25, eao, tk);
    
    # return effective Michaelis-Menten coefficient for CO2 
    return ( Kc * (1.0 + p->oi / Ko) ) ;
    
}