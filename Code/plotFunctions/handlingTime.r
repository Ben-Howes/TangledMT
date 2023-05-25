#######################################
## Plot handling time changes with mass
#######################################

library(tidyverse)

gpath = "/home/ben/Documents/TangledMT/Paper/Figures/"
setwd(gpath)

H0 = 1
B = 0.75
Rp = 0.5
a = 1
c = 0.5

masses = 1:100
masses = expand.grid(masses, masses) %>%
    rename("mi" = 1, "mj" = 2)

calculateHandling = function(mi, mj) {

    Hij = H0*(mj^(1-B))*(a*exp(-(((mj/mi) - Rp)^2)/2*(c)^2))

    return(Hij)

}

masses = masses %>% mutate(Hij = calculateHandling(mi, mj))

ggplot(filter(masses, mj == 100), aes(mi, Hij)) + geom_line()

ggplot(masses, aes(mj,mi, fill = Hij)) +
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