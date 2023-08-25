##############################################
## Create new formula for pR
## Where pR = pZ (probability of death)
## when H = 0
## So the growth rate is 0 at H = 0
##############################################

library(tidyverse)

gpath = "/home/ben/Documents/TangledMT/Paper/Figures/"
setwd(gpath)

x = seq(-100, 1000, length.out = 1001)
pZ = 0.2
masses = seq(log10(0.01), log10(100), length.out = 3)
masses = 10^masses
dat = expand.grid(x, masses) %>% 
    rename("x" = 1, "m" = 2)

pR = function(x, m, k = 0.1, pZ) {
    x0 = log((1/pZ)-1)/k
    G = 5*m^0.25
    pR = 1 / (1 + (exp(-k*(x - x0))))
    pR = pR/G
    return(pR)
}

dat = dat %>% mutate(G = m^0.25, pR = pR(x, m, k = 0.01, pZ = 0.2))

pRPlot = ggplot(dat, aes(x, pR, col = as.factor(round(log10(m), 2)))) + 
    geom_line(linewidth = 2.5) +
    theme_classic() +
    labs(x = expression(paste(H[i], "'")), y = expression(paste(p[R],sep="")),
    col = "Log10(Body\nMass)") +
    theme(text = element_text(size = 30)) +
    scale_colour_viridis_d()

pRPlot

ggsave(filename = "pRPlot.pdf", plot = pRPlot,  width = 18, height = 10)
