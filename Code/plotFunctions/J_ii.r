################################
## Plot interference term
################################

library(tidyverse)
library(lemon) ## facet_rep_wrap

gpath = "/home/ben/Documents/TangledMT/Paper/Figures/"
setwd(gpath)

H0 = 1
B = 0.75
k = 8.6173*(10^-5) ## Boltzmann constant
T0 = 293.15 ## 0 celsius in Kelvin
E = 0 ## Activation energy
V0 = 0.33
D0 = 1.62

temp = 20 ## teperatures in celsius
masses = seq(log10(0.01), log10(100), length.out = 3)
masses = 10^masses
N = seq(log10(0.000001), log10(1000000), length.out = 100)
N = 10^N
masses = expand.grid(masses, N, temp) %>% rename("mi" = 1, "Ni" = 2, "T" = 3)

calculateSearchRate = function(mi, mj, T) {

    aij = 2*(V0)*(D0)*(mi^(0.63))*(exp(1)^(-E/(k*(T + T0))))

    return(aij)

}

masses = masses %>% mutate(Jii = calculateSearchRate(mi, mi, 0), NiJii = Ni*mi*Jii)

NiJii = ggplot(masses, aes(log10(Ni), log10(NiJii), col = as.factor(log10(mi)), group = as.factor(log10(mi)))) +
    theme_classic() +
    geom_line(linewidth = 2.5) +
    scale_colour_viridis_d() +
    labs(x = "Log10(Consumer Density)", y = "Log10(Interference Rate Mass/Time)", col = "Log10\n(Consumer Mass)") +
    theme(text = element_text(size = 30))

NiJii

ggsave(filename = "NiJii.pdf", plot = NiJii,  width = 18, height = 10)
