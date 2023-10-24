#################################################
## Look at difference in mass distribution estimates
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
lnormFit = function(x) {

    mod = tryCatch(fitsad(x$M, "lnorm"),
    error = function(e) {NA})

    meanlog = tryCatch(coef(mod)[[1]], error = function(e) {NA})
    sdlog = tryCatch(coef(mod)[[2]], error = function(e) {NA})
    CV = tryCatch(coef(mod)[[2]]/coef(mod)[[1]], error = function(e) {NA})

    out = data.frame(seed = x$seed[[1]], g = x$g[[1]], meanlog, sdlog, CV)

    return(out)

}

dist = totalPopSpec %>%
    filter(pp == 0) %>%
    group_by(seed, g) %>%
    group_map(~lnormFit(.), .keep = TRUE) %>%
    bind_rows()

## Plot mean and sd of lognormal distribution over time
ggplot() +
    geom_rect(aes(xmin = 0, xmax = tClimate, ymin = -Inf, ymax = Inf), fill = "lightblue", alpha = 0.5) +
    geom_rect(aes(xmin = tClimate, xmax = tMax, ymin = -Inf, ymax = Inf), fill = "red", alpha = 0.25) +
    stat_lineribbon(data = dist, aes(g, meanlog), col = "black") +
    stat_lineribbon(data = dist, aes(g, sdlog), col = "grey50") +
    geom_vline(xintercept = tClimate, linetype = "dashed", linewidth = 2) +
    theme_classic() +
    labs(x = "Time Step", y = "Value") +
    theme(text = element_text(size = 30)) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_continuous(expand = c(0, 0))

ggsave(paste0("/home/ben/Documents/TangledMT/Paper/Figures/climateChange/climateMD.pdf"), width = 15, height = 10)
ggsave(paste0("/home/ben/Documents/TangledMT/Paper/Figures/climateChange/climateMD.png"), width = 15, height = 10)

## Look at difference in mean and sd of the log-normal distribution
## before and after climate change

compareLnorm = dist %>%
    mutate(clim = ifelse(g < 2.5e+08, 0, 1)) %>%
    filter(g > 2.0e+08 & g < 5.0e+08)

compareLnorm %>% 
    group_by(clim) %>%
    summarise(meanmeanlog = mean(meanlog, na.rm = T), meansdlog = mean(sdlog, na.rm = T), 
        varmeanlog = var(meanlog, na.rm = T), varsdlog = var(sdlog, na.rm = T))

## Plot different distributions before and after climate
avgFit = dist %>%
    filter(g >= 25000000) %>%
    mutate(climate = ifelse(g > tClimate, "Before Climate Change", "After Climate Change")) %>%
    group_by(climate) %>%
    summarise(meanlog = mean(meanlog), sdlog = mean(sdlog), CV = mean(CV))

## Make distributions before and after climate
climateDist = data.frame(beforeClimate = rlnorm(1000000, meanlog = avgFit$meanlog[1], sdlog = avgFit$sdlog[1]),
    afterClimate = rlnorm(1000000, meanlog = avgFit$meanlog[2], sdlog = avgFit$sdlog[2])) %>%
    pivot_longer(cols = c(beforeClimate, afterClimate), names_to = "climate", values_to = "value")

ggplot() +
    stat_halfeye(data = climateDist, aes(x = log10(value),y = climate, fill = climate, point_fill = climate),
    slab_colour = "black", slab_linewidth = 2, point_size = 7.5, shape = 21, stroke = 2.5, point_colour = "black",
    linewidth = 5) +
    theme_classic() +
    theme(legend.position = "none") +
    labs(x = "Log10(Mass)", y = NULL, fill = NULL, point_fill = NULL) +
    theme(text = element_text(size = 30)) +
    scale_y_discrete(labels = c("After\nClimage Change", "Before\nClimate Change"))

ggsave(paste0("/home/ben/Documents/TangledMT/Paper/Figures/climateChange/climateMDDist.pdf"), width = 15, height = 10)
ggsave(paste0("/home/ben/Documents/TangledMT/Paper/Figures/climateChange/climateMDDist.png"), width = 15, height = 10)