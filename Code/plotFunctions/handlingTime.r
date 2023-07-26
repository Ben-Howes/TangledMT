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
P0 = 1 ## Temperature constant (vary by taxa)
k = 8.6173*(10^-5) ## Boltzmann constant
T0 = 293.15 ## 0 celsius in Kelvin
E = 0 ## Activation energy

temp = 20 ## teperatures in celsius
masses = seq(log10(0.01), log10(100), length.out = 100)
masses = 10^masses
masses = expand.grid(masses, masses, temp) %>% rename("mi" = 1, "mj" = 2, "T" = 3)

calculateHandling = function(mi, mj, T) {

    Hij = H0*(mi^(-B))*(1-(a*exp(-(((mj/mi) - Rp)^2)/2*(c)^2)))*(exp(1)^(-E/(k*(T + T0))))

    return(Hij)

}

masses = masses %>% mutate(Hij = calculateHandling(mi, mj, T))

ggplot(masses, aes(log10(mi),log10(mj), fill = log10(Hij))) +
    geom_raster() +
    theme_classic() +
    scale_fill_viridis_c() +
    labs(x = "Log10 (Mass of Consumer)", y = "Log10 (Mass of Resource)", fill = "Log10 (Handling Time)") +
    theme(text = element_text(size = 30)) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_continuous(expand = c(0, 0)) +
    theme(legend.position = "bottom",
    legend.key.size = unit(1, 'cm'), #change legend key size
    legend.key.height = unit(1, 'cm'), #change legend key height
    legend.key.width = unit(3, 'cm'),
    legend.title = element_text(vjust = 1))

ggsave(filename = "handlingTime.png", width = 10, height = 10)

## Plot figure with 3 different consumer masses
masses = masses %>% filter(mi %in% c(min(masses$mi), masses[nrow(masses)/2, ]$mj, max(masses$mi))) %>%
    mutate(label = case_when(mi == min(masses$mi) ~ "Small Consumer", mi == masses[nrow(masses)/2, ]$mj ~ "Medium Consumer",
    mi == max(masses$mi) ~ "Large Consumer"))

handlingTime = ggplot(masses, aes(log10(mj), log10(Hij), col = as.factor(round(log10(mi), 1)))) +
    geom_line(linewidth = 5) +
    theme_classic() +
    labs(x = "Log10(Resource Body Mass)", y = "Log10(Handling Time)", col = "Log10(Consumer Body Mass)") +
    theme(strip.background = element_blank(),
    strip.text.x = element_blank(),
    text = element_text(size = 30),
    legend.position = "bottom") +
    scale_x_continuous(n.breaks = 8) +
    scale_colour_viridis_d()

handlingTime

ggsave(filename = "indHandlingTime.pdf", plot = handlingTime,  width = 18, height = 10)
