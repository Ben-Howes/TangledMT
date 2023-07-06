########################################
## Plot consumption rate function
########################################

library(tidyverse)
library(lemon) ## facet_rep_wrap

gpath = "/home/ben/Documents/TangledMT/Paper/Figures/"
setwd(gpath)

H0 = 1
B = 0.75
Rp = 0.1
a = 1
c = 1
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
masses = expand.grid(masses, masses, N, N, temp) %>% rename("mi" = 1, "mj" = 2, "Ni" = 3, "Nj" = 4, "T" = 5)

calculateAttackProb = function(mi, mj) {

    Aij = (1/(1 + 0.25*(exp(1)^-(mi^0.33))))*(1/(1 + (log10(Rp*(mi/mj))^2)))^5

    return(Aij)

}

calculateSearchRate = function(mi, mj, T) {

    aij = 2*(V0)*(D0)*(mi^(0.63))*(exp(1)^(-E/(k*(T + T0))))

    return(aij)

}

calculateHandling = function(mi, mj, T) {

    Hij = H0*(mi^(-B))*(1-(a*exp(-(((mj/mi) - Rp)^2)/2*(c)^2)))*(exp(1)^(-E/(k*(T + T0))))

    return(Hij)

}

consumptionRate = function(mi, mj, T, N) {

    Jij = (calculateSearchRate(mi, mj, T)*calculateAttackProb(mi, mj)*mj)/(1 + (calculateSearchRate(mi, mj, T)*calculateAttackProb(mi, mj)*calculateHandling(mi, mj, T)*N))
    return(Jij)

}

masses = masses %>% mutate(Jij = consumptionRate(mi, mj, T, Nj))

## Plot J_ij, which is per capita search rate
JijPlot = ggplot(masses, aes(log10(Nj), log10(Jij), col = as.factor(log10(mj)))) +
    facet_rep_wrap(. ~ log10(mi)) +
    geom_line(linewidth = 5) + 
    theme_classic() +
    labs(x = "Log10(Resource Density)", y = "Log10(Per-Capita Search Rate (Area/Time))", 
    col = "Log10(Resource\nBody Mass)", title = "Log10(Consumer Body Mass)") +
    theme(text = element_text(size = 30)) +
    scale_colour_viridis_d()

JijPlot

ggsave(filename = "Jij.pdf", plot = JijPlot,  width = 18, height = 10)

## Plot Jij multipled by resource density
NCjJijPlot = ggplot(masses, aes(log10(Nj), log10(Jij*Nj), col = as.factor(log10(mj)))) +
    facet_rep_wrap(. ~ log10(mi)) +
    geom_line(linewidth = 5) + 
    theme_classic() +
    labs(x = "Log10(Resource Density)", y = "Log10(Consumption Rate (Mass/Time))", 
    col = "Log10(Resource\nBody Mass)", title = "Log10(Consumer Body Mass)") +
    theme(text = element_text(size = 30)) +
    scale_colour_viridis_d()

NCjJijPlot

ggsave(filename = "NCjJij.pdf", plot = NCjJijPlot,  width = 18, height = 10)

## Now let's look at when our focal individual is the one being consumed

masses = masses %>% mutate(Jki = consumptionRate(mj, mi, T, Ni))

## Plot Jki, so when i, our focal individual, is the resource being consumed
NckJkiPlot = ggplot(filter(masses, Nj %in% c(1e-06, 1e+06)), aes(log10(Ni), log10(Jki*Nj), col = as.factor(log10(mi)))) +
    facet_rep_grid(log10(Nj) ~ log10(mj)) +
    geom_line(linewidth = 3) + 
    theme_classic() +
    labs(x = "Log10(Resource Density)", y = "Log10(Consumption Rate (Mass/Time))", 
    col = "Log10(Resource\nBody Mass)", title = "Cols = Log10(Consumer Body Mass)\nRows = Log10(Consumer Density)") +
    theme(text = element_text(size = 30)) +
    scale_colour_viridis_d()

NckJkiPlot

ggsave(filename = "NckJki.pdf", plot = NckJkiPlot,  width = 18, height = 10)
