#################################################
## Look at difference in population size
## and richness with climate
#################################################

library(tidyverse)
library(fitdistrplus) ## fit distributions to data
library(sads) ## fit logseries dist
library(ggdist)
library(parallel) ## mclapply

seeds = seq(100, 9900, 100)

## Load in totalPopSpec for each seed

getTotalPopSpec = function(x) {

    gpath = paste0("/home/ben/Documents/TangledMT/Results/climateChange/Seed_", x, "/Results/")
    setwd(gpath)

    totalPopSpec = read_delim("totalPopSpec.txt", col_names = FALSE) %>%
            rename(g = 1, s = 2, n = 3)
    traits = read_delim("../traits.txt", col_names = FALSE) %>%
            rename(M = 1, pp = 2) %>%
            add_column(.before = "M", s = 1:nrow(.))

    totalPop = totalPopSpec %>%
        group_by(g) %>%
        summarise(n = sum(n))

    ## Join trait data with totalpopspec
    totalPopSpec = totalPopSpec %>% left_join(traits)

    totalPopSpec = filter(totalPopSpec, pp == 0) %>%
        dplyr::select(g, s, n, M, pp) %>%
        mutate(.before = 1, seed = x)

    return(totalPopSpec)

}

totalPopSpec = lapply(seeds, getTotalPopSpec) 
totalPopSpec = totalPopSpec %>% bind_rows()

parameters = read_delim("/home/ben/Documents/TangledMT/Results/climateChange/Seed_100/Parameters.txt", col_names = FALSE) %>%
    rename(name = 1, var = 2, value = 3)
tClimate = parameters %>% filter(var == "climateT") %>% pull(value) %>% as.numeric()
tMax = max(totalPopSpec$g)

totalPop = totalPopSpec %>% 
    group_by(seed, g) %>%
    summarise(n = sum(n))

totalRich = totalPopSpec %>% 
    group_by(seed, g) %>%
    summarise(n = n())

## Population over time
ggplot() +
    geom_rect(aes(xmin = 0, xmax = tClimate, ymin = -Inf, ymax = Inf), fill = "lightblue", alpha = 0.5) +
    geom_rect(aes(xmin = tClimate, xmax = tMax, ymin = -Inf, ymax = Inf), fill = "red", alpha = 0.25) +
    stat_lineribbon(data = totalPop, aes(g, n)) +
    geom_vline(xintercept = tClimate, linetype = "dashed", linewidth = 2) +
    labs(x = "Time Step", y = "Total Population") +
    theme_classic() +
    theme(text = element_text(size = 20)) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_continuous(expand = c(0, 0))

ggsave(paste0("/home/ben/Documents/TangledMT/Paper/Figures/climateChange/climatePop.pdf"), width = 15, height = 10)

## Look at difference in population before and after climate change
## Looking at the population at 2.0e08 to 2.5e08, and 2.5e08 to 5.0e-08

comparePop = totalPop %>%
    mutate(clim = ifelse(g < 2.5e+08, 0, 1)) %>%
    filter(g > 2.0e+08 & g < 5.0e+08)

comparePop %>% 
    group_by(clim) %>%
    summarise(mean = mean(n))

## Richness over time
ggplot() +
    geom_rect(aes(xmin = 0, xmax = tClimate, ymin = -Inf, ymax = Inf), fill = "lightblue", alpha = 0.5) +
    geom_rect(aes(xmin = tClimate, xmax = tMax, ymin = -Inf, ymax = Inf), fill = "red", alpha = 0.25) +
    stat_lineribbon(data = totalRich, aes(g, n)) +
    geom_vline(xintercept = tClimate, linetype = "dashed", linewidth = 2) +
    labs(x = "Time Step", y = "Species Richness") +
    theme_classic() +
    theme(text = element_text(size = 20)) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_continuous(expand = c(0, 0))

ggsave(paste0("/home/ben/Documents/TangledMT/Paper/Figures/climateChange/climateRich.pdf"), width = 15, height = 10)


## Look at difference in population before and after climate change
## Looking at the population at 2.0e08 to 2.5e08, and 2.5e08 to 5.0e-08

compareRich = totalRich %>%
    mutate(clim = ifelse(g < 2.5e+08, 0, 1)) %>%
    filter(g > 2.0e+08 & g < 5.0e+08)

compareRich %>% 
    group_by(clim) %>%
    summarise(mean = mean(n))
