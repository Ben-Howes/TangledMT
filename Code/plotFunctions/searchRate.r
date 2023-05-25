#######################################
## Plot search rate changes with mass
#######################################

library(tidyverse)

gpath = "/home/ben/Documents/TangledMT/Paper/Figures/"
setwd(gpath)

V0 = 0.33
Pv = 0.21
D0 = 1.62
Pd = 0.21

masses = 1:100
masses = expand.grid(masses, masses) %>%
    rename("mi" = 1, "mj" = 2)

calculateSearchRate = function(mi, mj) {

    aij = 2*(V0)*(D0)*(mj^(Pv + Pd - 1))*(mi^(Pd))

    return(aij)

}

masses = masses %>% mutate(aij = calculateSearchRate(mi, mj))

ggplot(masses, aes(mj,mi, fill = aij)) +
    geom_raster() +
    theme_classic() +
    scale_fill_viridis_c() +
    labs(x = "Mass of Consumer", y = "Mass of Resource", fill = "Search Rate") +
    theme(text = element_text(size = 30)) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_continuous(expand = c(0, 0)) +
    theme(legend.position = "bottom",
    legend.key.size = unit(1, 'cm'), #change legend key size
    legend.key.height = unit(1, 'cm'), #change legend key height
    legend.key.width = unit(3, 'cm'),
    legend.title = element_text(vjust = 1))

ggsave(filename = "searchRate.png", width = 10, height = 10)
