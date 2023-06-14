######################################
## Plot the logistic growth rate
## used for primary producers
######################################

library(tidyverse)

Mi = 0.1 ## Mass of focal primary producer i
r = Mi^(-0.15) ## Intrinsic growth rate of primary producer i
Mc = seq(log10(0.001), log10(100), length.out = 100) ## Total mass of primary producers in cell C
Mc = Mi + (10^Mc) ## Transform from log10, and add mass of Mi, since the total mass cannot be less than this
K = 100 ## Constant carry capacity

logGrowth = function(r, M) {

    growth = r*(1-(Mc/K))
    return(growth)

}

results = data.frame(Mc) %>% mutate(g = logGrowth(r, Mc))

ggplot(results, aes(log10(Mc), g)) + 
    geom_point() + 
    theme_classic()
