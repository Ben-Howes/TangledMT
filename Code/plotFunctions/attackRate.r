#######################################
## Plot attack probability changes with mass
#######################################

library(tidyverse)

gpath = "/home/ben/Documents/TangledMT/Paper/Figures/"
setwd(gpath)

Rp = 0.1

masses = 1:100
masses = expand.grid(masses, masses) %>%
    rename("mi" = 1, "mj" = 2)

calculateAttackProb = function(mi, mj) {

    Aij = (1/(1 + 0.25*(exp(1)^-(mj^0.33))))*(1/(1 + (log10(Rp*(mj/mi))^2)))^0.2

    return(Aij)

}

masses = masses %>% mutate(Aij = calculateAttackProb(mi, mj))

ggplot(masses, aes(mj,mi, fill = Aij)) +
    geom_raster() +
    theme_classic() +
    scale_fill_viridis_c() +
    labs(x = "Mass of Consumer", y = "Mass of Resource", fill = "Attack Probability") +
    theme(text = element_text(size = 30)) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_continuous(expand = c(0, 0)) +
    theme(legend.position = "bottom",
    legend.key.size = unit(1, 'cm'), #change legend key size
    legend.key.height = unit(1, 'cm'), #change legend key height
    legend.key.width = unit(3, 'cm'),
    legend.title = element_text(vjust = 1))

ggsave(filename = "attackProbability.png", width = 10, height = 10)
