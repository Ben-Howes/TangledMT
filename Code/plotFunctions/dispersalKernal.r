#####################################
## Make example dispersal kernals
## given mass
#####################################

library(tidyverse)
library(lemon) ## facet_rep_wrap

gpath = "/home/ben/Documents/TangledMT/Paper/Figures/"
setwd(gpath)

masses = data.frame(Mi = seq(log10(0.01), log10(100), length.out = 3))
masses = masses %>% mutate(Mi = 10^Mi)

## Now we need to make a distribution of dispersals
## probabilities given this data
## using the negative exponential distribution as suggested by sutherland

probDisp = function(x, Mi) {
    probDisp = exp(1)^(-x/1.5)
    return(probDisp)
}

density = expand.grid(x = rexp(100000, 1), masses$Mi) %>%
    rename(x = 1, Mi = 2) %>%
    mutate(x = x*(Mi^0.63))

ggplot(density, aes(x, col = as.factor(log10(Mi)))) +
    facet_rep_wrap(. ~ Mi, scales = "free_x") +
    geom_density(linewidth = 2) +
    theme_classic() +
    theme(text = element_text(size = 30)) +
    labs(x = "Distance (D)", y = "Density", col = "Log10\n(Body Mass)") +
    scale_colour_viridis_d()

ggsave(filename = "dispersalDensity.png", width = 15, height = 10)

data = expand.grid(seq(0, 100, 0.01), masses$Mi) %>%
    rename(x = 1, Mi = 2) %>%
    mutate(x = x*(Mi^0.63), pDisp = probDisp(x, Mi))

ggplot(data, aes(x, pDisp, col = as.factor(log10(Mi)))) +
    facet_rep_wrap(. ~ Mi, scales = "free_x") +
    geom_line(linewidth = 2) +
    theme_classic() + 
    labs(x = "Distance (D)", y = "Probability of Dispersal > D", col = "Log10\n(Body Mass)") +
    theme(text = element_text(size = 30)) +
    scale_colour_viridis_d()

ggsave(filename = "dispersalKernal.png", width = 15, height = 10)

testFunc = function(x) {
    out = exp(1)^(-1.5*x)
    return(out)
}

test = data.frame(x = seq(0, 100, by = 0.01)) %>%
    mutate(y = testFunc(x))

ggplot(test, aes(x, y)) + 
    geom_line()
