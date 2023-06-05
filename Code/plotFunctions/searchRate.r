#######################################
## Plot search rate changes with mass
#######################################

library(tidyverse)

gpath = "/home/ben/Documents/TangledMT/Paper/Figures/"
setwd(gpath)

V0 = 0.33
D0 = 1.62
P0 = 1 ## Temperature constant (vary by taxa)
k = 8.6173*(10^-5) ## Boltzmann constant
T0 = 293.15 ## 0 celsius in Kelvin
E = 0 ## Activation energy

temp = 20 ## teperatures in celsius
masses = seq(log10(0.0000001), log10(1000000), length.out = 100)
masses = 10^masses
masses = expand.grid(masses, masses, temp) %>% rename("mi" = 1, "mj" = 2, "T" = 3)

calculateSearchRate = function(mi, mj, T) {

    aij = 2*(V0)*(D0)*(mi^(0.63))*(mj^(0.21))*(P0)*(exp(1)^(-E/(k*(T + T0))))

    return(aij)

}

masses = masses %>% mutate(aij = calculateSearchRate(mi, mj, T))

ggplot(masses, aes(log10(mi), log10(mj), fill = log10(aij))) +
    geom_raster() +
    theme_classic() +
    scale_fill_viridis_c() +
    labs(x = "Log10 (Body Mass of Consumer)", y = "Log10 (Body Mass of Resource)", fill = "Log10 (Search Rate)") +
    theme(text = element_text(size = 30)) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_continuous(expand = c(0, 0)) +
    theme(legend.position = "bottom",
    legend.key.size = unit(1, 'cm'), #change legend key size
    legend.key.height = unit(1, 'cm'), #change legend key height
    legend.key.width = unit(3, 'cm'),
    legend.title = element_text(vjust = 1))

ggsave(filename = "searchRate.png", width = 10, height = 10)
