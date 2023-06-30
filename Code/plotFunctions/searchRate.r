#######################################
## Plot search rate changes with mass
#######################################

library(tidyverse)
library(lemon) ## facet_rep_wrap

gpath = "/home/ben/Documents/TangledMT/Paper/Figures/"
setwd(gpath)

V0 = 0.33
D0 = 1.62
P0 = 1 ## Temperature constant (vary by taxa)
k = 8.6173*(10^-5) ## Boltzmann constant
T0 = 293.15 ## 0 celsius in Kelvin
E = 0 ## Activation energy
Rp = 0.1

temp = 20 ## teperatures in celsius
masses = seq(log10(0.000001), log10(1000000), length.out = 100)
masses = 10^masses
masses = expand.grid(masses, masses) %>% rename("mi" = 1, "mj" = 2)

calculateAttackProb = function(mi, mj) {

    Aij = (1/(1 + 0.25*(exp(1)^-(mi^0.33))))*(1/(1 + (log10(Rp*(mi/mj))^2)))^0.2

    return(Aij)

}

calculateSearchRate = function(mi, mj, T) {

    aij = 2*(V0)*(D0)*(mi^(0.63))*calculateAttackProb(mi, mj)*(exp(1)^(-E/(k*(T + T0))))

    return(aij)

}

masses = masses %>% mutate(aij = calculateSearchRate(mi, mj, T))

## Plot figure with 3 different consumer masses
masses = masses %>% filter(mi %in% c(min(masses$mi), masses[nrow(masses)/2, ]$mj, max(masses$mi))) %>%
    mutate(label = case_when(mi == min(masses$mi) ~ "Small Consumer", mi == masses[nrow(masses)/2, ]$mj ~ "Medium Consumer",
    mi == max(masses$mi) ~ "Large Consumer"))

searchRate = ggplot(masses, aes(log10(mj), log10(aij), col = as.factor(round(log10(mi), 1)))) +
    geom_line(linewidth = 5) +
    theme_classic() +
    labs(x = "Log10(Resource Body Mass)", y = "Log10(Search Rate)", col = "Log10(Consumer Body Mass)") +
    theme(strip.background = element_blank(),
    strip.text.x = element_blank(),
    text = element_text(size = 30),
    legend.position = "bottom") +
    scale_x_continuous(n.breaks = 8) +
    scale_colour_viridis_d()

searchRate

ggsave(filename = "searchRate.png", plot = searchRate,  width = 18, height = 10)
