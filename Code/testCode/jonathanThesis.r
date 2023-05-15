####################################################
## First try at creating interaction strengths
## from MTE based on equations in
## Jonathan thesis
####################################################

library(tidyverse)

## Constant set from thesis
b0 = 1.71e-6 ## Normalising constant for birth rate
d0 = 4.15e-8 ## Normalising constant for death rate
B = 0.75 ## Scaling exponenet for metabolism
a0 = 8.31e-4 ## Normalising constant for search rate
pv = 0.26 ## Scaling exponent for velocity
pd = 0.21 ## scaling exponent for detection distance
# e01 = 0.2 ## biomass conversion efficient, feeding on plant
# e02 = 0.5 ## biomass conversion efficient, feeding on animal

## mass of theoretical species (kg)
m1 = 1
m2 = 10

interactingScaling = function(i, j, captureType = c("active", "grazing")){

    if(captureType == "active") {
        int = (j^(pv + (2*pd)))*(sqrt(1 + ((i/j)^(2*pv))*((i/j)^pd)))
    }

    if(captureType == "grazing") {
        int = (j*(pv + (2*pd)))*((i/j)^pd)
    }

    return(int)
}


interactStrength = function(i, j, captureType, interactType = c("prey", "predator")) {

    if(interactType == "prey") {
        int = -a0*interactingScaling(i,j, captureType = captureType)*(j^-1)
    }

    if(interactType == "predator") {
        if(captureType == "grazing") {e0 = 0.2} else {e0 = 0.5}
        int = e0*a0*interactingScaling(j,i, captureType = captureType)*(i^-1)
    }

return(int)

}

interactStrength(m1, m2, captureType = "active", interactType = "predator")


## Now let's make lots of different masses and see what values for interaction strength we get

masses = seq(0.5, 10, by = 0.5)
masses = expand.grid(masses, masses) %>% rename("m1" = 1, "m2" = 2)

test = masses %>% mutate(interact = interactStrength(m1, m2, captureType = "active", interactType = "predator"))

ggplot(test, aes(m1, interact, col = as.factor(m2))) + geom_line(linewidth = 2) + theme_classic()
