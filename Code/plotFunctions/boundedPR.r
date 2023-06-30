##############################################
## Create new formula for Poff
## Where Poff = 0 at H < 0 and 1 at H > 1
##############################################

library(tidyverse)

gpath = "/home/ben/Documents/TangledMT/Paper/Figures/"
setwd(gpath)

dat = data.frame(x = seq(0, 1, length.out = 101))

pOff = function(x, k = 10) {
    pOff = 1 / (1 + exp(-k * (x - 0.5)))
    return(pOff)
}

dat = dat %>% mutate(pOff = pOff(x))

pOffPlot = ggplot(dat, aes(x, pOff)) + 
    geom_line(linewidth = 2.5) +
    theme_classic() +
    labs(x = expression(paste(H[i], "'")), y = expression(paste(p[R],sep=""))) +
    theme(text = element_text(size = 30))

pOffPlot

ggsave(filename = "pOffPlot.png", plot = pOffPlot,  width = 18, height = 10)
