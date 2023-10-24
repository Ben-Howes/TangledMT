#################################################
## Look at difference in SAD estimates
## with climate change
#################################################

library(tidyverse)
library(ggdist)
library(parallel) ## mclapply
library(sads) ## fit logseries dist

seeds = seq(100, 9900, 100)

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

## Plot log-series estimates over time
lsFit = function(x) {

    mod = tryCatch(fitsad(x$n, "ls"),
    error = function(e) {NA})

    N = tryCatch(coef(mod)[[1]], error = function(e) {NA})
    alpha = tryCatch(coef(mod)[[2]], error = function(e) {NA})
    alphaN = tryCatch(coef(mod)[[1]]/coef(mod)[[2]], error = function(e) {NA})

    out = data.frame(seed = x$seed[[1]], g = x$g[[1]], alpha, N, alphaN)

    return(out)

}

fit = totalPopSpec %>% 
    group_split(seed, g) %>%
    mclapply(., lsFit, mc.cores = 8) %>% 
    bind_rows()

## Alpha
ggplot() +
    geom_rect(aes(xmin = 0, xmax = tClimate, ymin = -Inf, ymax = Inf), fill = "lightblue", alpha = 0.5) +
    geom_rect(aes(xmin = tClimate, xmax = tMax, ymin = -Inf, ymax = Inf), fill = "red", alpha = 0.25) +
    stat_lineribbon(data = fit, aes(g, alpha)) +
    geom_vline(xintercept = tClimate, linetype = "dashed", linewidth = 2) +
    labs(x = "Time Step", y = "Log-Series Alpha") +
    theme_classic() +
    theme(text = element_text(size = 20)) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_continuous(expand = c(0, 0))

ggsave(paste0("/home/ben/Documents/TangledMT/Paper/Figures/climateChange/climateAlpha.pdf"), width = 15, height = 10)

## Look at difference in alpha before and after climate change
## Looking at the population at 2.0e08 to 2.5e08, and 2.5e08 to 5.0e-08

compareAlpha = fit %>%
    mutate(clim = ifelse(g < 2.5e+08, 0, 1)) %>%
    filter(g > 2.0e+08 & g < 5.0e+08)

compareAlpha %>% 
    group_by(clim) %>%
    summarise(mean = mean(alpha, na.rm = T))


ggplot() +
    geom_rect(aes(xmin = 0, xmax = tClimate, ymin = -Inf, ymax = Inf), fill = "lightblue", alpha = 0.5) +
    geom_rect(aes(xmin = tClimate, xmax = tMax, ymin = -Inf, ymax = Inf), fill = "red", alpha = 0.25) +
    stat_lineribbon(data = fit, aes(g, alphaN)) +
    geom_vline(xintercept = tClimate, linetype = "dashed", linewidth = 2) +
    labs(x = "Time Step", y = "Log-Series Alpha/N") +
    theme_classic() +
    theme(text = element_text(size = 20)) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_continuous(expand = c(0, 0))

ggsave(paste0("/home/ben/Documents/TangledMT/Paper/Figures/climateChange/climateAlphaN.pdf"), width = 15, height = 10)