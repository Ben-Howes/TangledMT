#####################################
## Make example dispersal kernals
## given mass
#####################################

gpath = "/home/ben/Documents/TangledMT/Paper/Figures/"
setwd(gpath)

masses = data.frame(Mi = seq(log10(0.01), log10(100), length.out = 3))
masses = masses %>% mutate(Mi = 10^Mi)

## Find normalisation constant so smallest species always has maximum dispersal of 1
D0 = 1/(min(masses)^0.63)

## Maximum dispersal functiion without temperature
maxDispersal = function(Mi) {
    maxDisp = D0*Mi^0.63
}

## Calculate maximum dispersal for each species
masses = masses %>% mutate(disp = maxDispersal(Mi))

## Now we need to make a distribution of dispersals
## probabilities given this data
## From reading, the Weibull distribution is best

x = data.frame(x = rweibull(10000, 10))
ggplot(x, aes(x)) + geom_density()
