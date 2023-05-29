#######################################
## Plot handling time changes with mass
#######################################

library(tidyverse)

gpath = "/home/ben/Documents/TangledMT/Paper/Figures/"
setwd(gpath)

H0 = 1
B = 0.75
Rp = 0.1
a = 1
c = 1

masses = seq(log10(0.00001), log10(1000000), length.out = 100)
masses = 10^masses
masses = expand.grid(masses, masses) %>%
    rename("mi" = 1, "mj" = 2)

calculateHandling = function(mi, mj) {

    Hij = H0*(mi^(B))*(1-(a*exp(-(((mj/mi) - Rp)^2)/2*(c)^2)))

    return(Hij)

}

masses = masses %>% mutate(Hij = calculateHandling(mi, mj))

ggplot(filter(masses, mi == max(masses$mi)), aes(mj, Hij)) + geom_line()

ggplot(masses, aes(log10(mi),log10(mj), fill = log10(Hij))) +
    geom_raster() +
    theme_classic() +
    scale_fill_viridis_c() +
    labs(x = "Mass of Consumer", y = "Mass of Resource", fill = "Handling Time") +
    theme(text = element_text(size = 30)) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_continuous(expand = c(0, 0)) +
    theme(legend.position = "bottom",
    legend.key.size = unit(1, 'cm'), #change legend key size
    legend.key.height = unit(1, 'cm'), #change legend key height
    legend.key.width = unit(3, 'cm'),
    legend.title = element_text(vjust = 1))

ggsave(filename = "handlingTime.png", width = 10, height = 10)