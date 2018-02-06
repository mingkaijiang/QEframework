assim <- function(ci, gamma_star, a1, a2) {
    #    Morning and afternoon calcultion of photosynthesis with the
    #    limitation defined by the variables passed as a1 and a2, i.e. if we
    #     are calculating vcmax or jmax limited.
    
    #    Parameters:
    #    ----------
    #    ci : float
    #    intercellular CO2 concentration.
    #    gamma_star : float
    #    CO2 compensation point in the abscence of mitochondrial respiration
    #    a1 : float
    #    variable depends on whether the calculation is light or rubisco
    #    limited.
    #    a2 : float
    #    variable depends on whether the calculation is light or rubisco
    #    limited.
    
    #    Returns:
    #    -------
    #    assimilation_rate : float
    #    assimilation rate assuming either light or rubisco limitation.

        if (ci < gamma_star)
            return (0.0)
    else
        return (a1 * (ci - gamma_star) / (a2 + ci));
    
}