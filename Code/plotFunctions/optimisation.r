#######################################
## Create optimal foragning function
## and plot example output
#######################################


library(tidyverse)

gpath = "/home/ben/Documents/TangledMT/Paper/Figures/"
setwd(gpath)

#####
## Set all constants needed for the functions
#####

V0 = 0.33
D0 = 1.62
P0 = 1 ## Temperature constant (vary by taxa)
k = 8.6173*(10^-5) ## Boltzmann constant
T0 = 293.15 ## 0 celsius in Kelvin
E = 0 ## Activation energy
H0 = 1
B = 0.75
Rp = 0.1
a = 1
c = 1

#####
## Set mass and temperature for some individuals
#####

temp = 20 ## teperatures in celsius
masses = seq(log10(0.0000001), log10(1000000), length.out = 1000)
masses = 10^masses
masses = expand.grid(masses, masses, temp) %>% rename("mi" = 1, "mj" = 2, "T" = 3)

## Write in functions needed

calculateSearchRate = function(mi, mj, T) {

    aij = 2*(V0)*(D0)*(mi^(0.63))*(mj^(0.21))*(P0)*(exp(1)^(-E/(k*(T + T0))))

    return(aij)

}

calculateHandling = function(mi, mj, T) {

    Hij = H0*(mi^(-B))*(1-(a*exp(-(((mj/mi) - Rp)^2)/2*(c)^2)))*(P0)*(exp(1)^(-E/(k*(T + T0))))

    return(Hij)

}


calculateAttackProb = function(mi, mj) {

    Aij = (1/(1 + 0.25*(exp(1)^-(mi^0.33))))*(1/(1 + (log10(Rp*(mi/mj))^2)))^0.2

    return(Aij)

}

calculateProfitability = function(mi, mj, T) {

    pij = (calculateAttackProb(mi, mj)*mj)/calculateHandling(mi, mj, T)

    return(pij)

}

#####
## Choose a focal consumer species
#####

masses = masses %>% filter(mi == masses[nrow(masses)/2,]$mj)

## Find most profitable species for our focal species
## and arrange from most to least profitable

masses = masses %>% mutate(pij = calculateProfitability(mi, mj, T)) %>%
    arrange(-pij)

## Add individuals assuming they are all equally abundant
masses = masses %>% mutate(N = mj^-0.75)

optimisation = function(dat, x) {

    mi = dat[1,]$mi
    mj = dat[x,]$mj

    upper = sum(calculateSearchRate(mi, dat[1:x,]$mj, T)*calculateAttackProb(mi, dat[1:x,]$mj)*dat[1:x,]$N*dat[1:x,]$mj)
    lower = 1 + sum(calculateSearchRate(mi, dat[1:x,]$mj, dat[1:x, ]$T)*calculateAttackProb(mi, dat[1:x,]$mj)*calculateHandling(mi, dat[1:x,]$mj, dat[1:x,]$T)*dat[1:x,]$N)

    out = data.frame(mj, pij = dat[x,]$pij, opt = upper/lower)

    return(out)

}

opt = lapply(1:nrow(masses), function(x) optimisation(x = x, dat = masses)) %>% bind_rows()

ggplot(opt, aes(x = pij, opt)) + 
    geom_point(size = 7.5, aes(col = log10(mj))) +
    theme_classic() +
    labs(x = "Profitability", y = "Consumption Rate (kg/s)", 
    col = "Log10\n(Resouce Body Mass)") +
    scale_colour_viridis_c() +
    theme(text = element_text(size = 30),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank()) +
    theme(legend.position = "bottom",
    legend.key.size = unit(1, 'cm'), #change legend key size
    legend.key.height = unit(1, 'cm'), #change legend key height
    legend.key.width = unit(3, 'cm'),
    legend.title = element_text(vjust = 1)) +
    scale_x_reverse()

ggsave(filename = "optimisation.png", width = 10, height = 10)
