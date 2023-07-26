##############################################
## Create new formula for Poff
## Where Poff = 0 at H < 0 and 1 at H > 1
##############################################

library(tidyverse)

gpath = "/home/ben/Documents/TangledMT/Paper/Figures/"
setwd(gpath)

x = seq(-10, 10, length.out = 101)
masses = seq(log10(0.01), log10(100), length.out = 3)
masses = 10^masses
dat = expand.grid(x, masses) %>% rename("x" = 1, "m" = 2)

pOff = function(x, m, k = 0.1) {
    G = 5*m^0.25
    pOff = 1 / (1 + (exp(-(x - 0.5))))
    pOff = pOff/G
    return(pOff)
}

dat = dat %>% mutate(G = m^0.25, pOff = pOff(x, m, k = 1))

pOffPlot = ggplot(dat, aes(x, pOff, col = as.factor(round(log10(m), 2)))) + 
    geom_line(linewidth = 2.5) +
    theme_classic() +
    labs(x = expression(paste(H[i], "'")), y = expression(paste(p[R],sep="")),
    col = "Log10(Body\nMass)") +
    theme(text = element_text(size = 30)) +
    scale_colour_viridis_d()

pOffPlot

ggsave(filename = "pOffPlot.pdf", plot = pOffPlot,  width = 18, height = 10)
