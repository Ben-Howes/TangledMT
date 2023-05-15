###########################################
## Create functions for ADBM from
## Petchley 2008
###########################################

library(tidyverse)
library(viridis)

set.seed(123)

## Mi is the resource
## Mj is the consumer

## All constants used are taken from the MS where they use different parameters per food web
## I set this data to the Coachella dataset which is one of the terrestrial datasets whcih was predicted best by the ADBM
## e, h, and n are all arbitrarily set in the MS to 1, and here as well

###########################################
## Handling time function
###########################################
## Handling time will always be Inf if (Mi/Mj) >= b

handlingTime = function(Mi, Mj, b = 0.6, h = 1) {
    if((Mi/Mj) < b) {
        H = (h/(b - (Mi/Mj)))
    } else {
        H = NA
    }

    return(H)
}

## Vectorize for use in mutate
handlingTime = Vectorize(handlingTime)

## Create data frame of masses to see relationship
masses = seq(1,10,by=0.1)
masses = expand.grid(masses, masses) %>% rename("Mi" = 1, "Mj" = 2)

allDat = masses %>% mutate(H = handlingTime(Mi, Mj))

ggplot(allDat, aes(Mj, H, col = as.factor(Mi))) + 
    geom_line(linewidth = 2) +
    theme_classic() +
    labs(x = "Consumer Mass", y = "Handling Time", col = "Resource Mass") +
    theme(text = element_text(size = 20)) +
        scale_colour_viridis_d()

###########################################
## Resource Density Function
###########################################

## Constants
n = 1
ni = -0.75

resourceDensity = function(Mi) {
    N = n*(Mi^ni)
    return(N)
}

allDat = allDat %>% mutate(N = resourceDensity(Mi))

ggplot(data.frame(Mi = seq(1, 10, by = 0.1), N = resourceDensity(seq(1,10,by = 0.1))), aes(Mi, N)) +
    geom_line(linewidth = 2) +
    theme_classic() +
    labs(x = "Resource Mass", y = "Resource Density") +
    theme(text = element_text(size = 20))  

###########################################
## Energy Content Function
###########################################

## Constants
## In Jonathan Thesis he uses e = 0.2 for grazing, and e = 0.5 for consumers
e = 1

energyContent = function(Mi) {
    E = e*(Mi)
    return(E)
}

allDat = allDat %>% mutate(E = energyContent(Mi))

ggplot(data.frame(Mi = seq(1, 10, by = 0.1), E = energyContent(seq(1,10,by = 0.1))), aes(Mi, E)) +
    geom_line(linewidth = 2) +
    theme_classic() +
    labs(x = "Resource Mass", y = "Energy Content") +
    theme(text = element_text(size = 20))  

###########################################
## Attack Rate Function
###########################################

## Constants
a = 1
ai = 0.58
aj = 0.33

attackRate = function(Mi, Mj) {
    A = a*(Mi^ai)*(Mj^aj)
    return(A)
}

allDat = allDat %>% mutate(A = attackRate(Mi, Mj))

ggplot(allDat, aes(Mj, A, col = as.factor(Mi))) +
    geom_line(linewidth = 2) +
    theme_classic() +
    labs(x = "Consumer Mass", y = "Attack Rate", col = "Resource Mass") +
    theme(text = element_text(size = 20)) +
    scale_colour_viridis_d()

###########################################
## Encounter Rate Function
###########################################

## Encounter rate is just the product of attack rate and density
encounterRate = function(Mi, Mj) {
    Y = resourceDensity(Mi)*attackRate(Mi, Mj)
    return(Y)
}

allDat = allDat %>% mutate(Y = A*N)

ggplot(allDat, aes(Mj, Y, col = as.factor(Mi))) +
    geom_line(linewidth = 2) +
    theme_classic() +
    labs(x = "Consumer Mass", y = "Encounter Rate", col = "Resource Mass") +
    theme(text = element_text(size = 20)) +
    scale_colour_viridis_d()

###########################################
## Profitability Function
###########################################

## For each consumer J we can calculate the profitability of feeding on a resource I

profitability = function(Mi, Mj) {
    P = energyContent(Mi)/handlingTime(Mi, Mj)
    return(P)
}

allDat = allDat %>% mutate(P = profitability(Mi, Mj))

ggplot(allDat, aes(Mj, P, col = as.factor(Mi))) + 
    geom_line(linewidth = 2) + 
    theme_classic() +
    labs(x = "Consumer Mass", y = "Profitability", col = "Resource Mass") +
    theme(text = element_text(size = 20)) +
    scale_colour_viridis_d()

###########################################
## Optimal Foraging Function
###########################################

## The final model combines all functions to predict the
## optimal diet k of each consumer J that maximises energy intake
## We start by ranking resources I by profitability for consumer J in decreasing order
## Then we add decreasingly less profitable resources I to our optimsation function
## Until our optimal foragning value begins to decrease

optimalForaging = function(dat) {

## Calculate Profitability
dat = dat %>% mutate(P = profitability(Mi, Mj))

## Rank by Profitability P
dat = dat %>% arrange(-P)

OF = 0
optimalForaging = NULL
for(x in 1:nrow(dat)) {
    OFnew = sum(encounterRate(dat$Mi[1:x], dat$Mj[1:x])*energyContent(dat$Mi[1:x]))/(1+sum(encounterRate(dat$Mi[1:x], dat$Mj[1:x])*handlingTime(dat$Mi[1:x], dat$Mj[1:x])))
    if(is.na(OFnew)){break} else if(OFnew > OF) {
        OF = OFnew
        optimalForaging = bind_rows(optimalForaging, 
        add_column(dat[x,], F = sum(encounterRate(dat$Mi[x], dat$Mj[x])*energyContent(dat$Mi[x]))/(1+sum(encounterRate(dat$Mi[x], dat$Mj[x])*handlingTime(dat$Mi[x], dat$Mj[x]))),
        OF))} else {
            break
            }
}

return(optimalForaging)

}

## Calculate optimal diet for all consumers in our allDat data frame
splitDat = allDat %>% group_split(Mj)
optimalForagingSets = lapply(splitDat, optimalForaging) %>% bind_rows()

ggplot(optimalForagingSets, aes(Mj, Mi, fill = F)) + 
    geom_raster() +
    theme_classic() +
    theme(text = element_text(size = 20)) +
    labs(x = "Consumer Mass", y = "Resource Mass", col = "Rate of\nEnergy Intake") +
    scale_fill_viridis_c()

#### Run for some larger networks with more variation in body mass
newMasses = rlnorm(15, meanlog = 1, sdlog = 5)
newMasses = expand.grid(newMasses, newMasses) %>% rename("Mi" = 1, "Mj" = 2)
splitNewMasses = newMasses %>% group_split(Mj)
testOF = lapply(splitNewMasses, optimalForaging) %>% bind_rows()

ggplot(testOF, aes(Mj, Mi, col = F)) + 
    geom_point(size = 10) +
    theme_classic() +
    theme(text = element_text(size = 20)) +
    labs(x = "Consumer Mass", y = "Resource Mass", col = "Rate of\nEnergy Intake") +
    scale_colour_viridis_c()

testWide = testOF %>% pivot_wider(id_cols = c(Mi), values_from = F, names_from = Mj, values_fill = 0)
