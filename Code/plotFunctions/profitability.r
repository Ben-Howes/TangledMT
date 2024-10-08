#######################################
## Plot profitability
######################################

library(tidyverse)

gpath = "/home/ben/Documents/TangledMT/Paper/Figures/"
setwd(gpath)

H0 = 1
B = 0.75
Rp = 0.1
a = 1
c = 1
P0 = 1 ## Temperature constant (vary by taxa)
k = 8.6173*(10^-5) ## Boltzmann constant
T0 = 293.15 ## 0 celsius in Kelvin
E = 0 ## Activation energy

temp = 20 ## teperatures in celsius
masses = seq(log10(0.0000001), log10(1000000), length.out = 100)
masses = 10^masses
masses = expand.grid(masses, masses, temp) %>% rename("mi" = 1, "mj" = 2, "T" = 3)

calculateAttackProb = function(mi, mj) {

    Aij = (1/(1 + 0.25*(exp(1)^-(mi^0.33))))*(1/(1 + (log10(Rp*(mi/mj))^2)))^0.2

    return(Aij)

}

calculateHandling = function(mi, mj, T) {

    Hij = H0*(mi^(-B))*(1-(a*exp(-(((mj/mi) - Rp)^2)/2*(c)^2)))*(P0)*(exp(1)^(-E/(k*(T + T0))))

    return(Hij)

}

calculateProfitability = function(mi, mj, T) {

    pij = (calculateAttackProb(mi, mj)*mj)/calculateHandling(mi, mj, T)

    return(pij)

}

masses = masses %>% mutate(pij = calculateProfitability(mi, mj, T))

ggplot(masses, aes(log10(mj), log10(mi), fill = log10(pij))) +
    geom_raster() +
    theme_classic() +
    scale_fill_viridis_c() +
    labs(x = "Log10 (Body Mass of Consumer)", y = "Log10 (Body Mass of Resource)", 
    fill = "Log10 (Profitability)") +
    theme(text = element_text(size = 30)) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_continuous(expand = c(0, 0)) +
    theme(legend.position = "bottom",
    legend.key.size = unit(1, 'cm'), #change legend key size
    legend.key.height = unit(1, 'cm'), #change legend key height
    legend.key.width = unit(3, 'cm'),
    legend.title = element_text(vjust = 1))

ggsave(filename = "profitability.png", width = 10, height = 10)
