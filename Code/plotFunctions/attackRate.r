#######################################
## Plot attack probability changes with mass
#######################################

library(tidyverse)

gpath = "/home/ben/Documents/TangledMT/Paper/Figures/"
setwd(gpath)

Rp = 0.1

masses = seq(log10(0.0000001), log10(1000000), length.out = 100)
masses = 10^masses
masses = expand.grid(masses, masses) %>% rename("mi" = 1, "mj" = 2)


calculateAttackProb = function(mi, mj) {

    Aij = (1/(1 + 0.25*(exp(1)^-(mi^0.33))))*(1/(1 + (log10(Rp*(mi/mj))^2)))^0.2

    return(Aij)

}

masses = masses %>% mutate(Aij = calculateAttackProb(mi, mj))

ggplot(masses, aes(log10(mj), log10(mi), fill = Aij)) +
    geom_raster() +
    theme_classic() +
    scale_fill_viridis_c() +
    labs(x = "Log10 (Body Mass of Consumer)", y = "Log10 (Body Mass of Resource)", fill = "Attack Probability") +
    theme(text = element_text(size = 30)) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_continuous(expand = c(0, 0)) +
    theme(legend.position = "bottom",
    legend.key.size = unit(1, 'cm'), #change legend key size
    legend.key.height = unit(1, 'cm'), #change legend key height
    legend.key.width = unit(3, 'cm'),
    legend.title = element_text(vjust = 1))

ggsave(filename = "attackProbability.png", width = 10, height = 10)
