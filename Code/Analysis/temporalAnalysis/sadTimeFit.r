
#################################################
## Calculate Alpha, N, and alphaN
## over time with mean and sd
#################################################

library(tidyverse)
library(fitdistrplus) ## fit distributions to data
library(sads) ## fit logseries dist
library(ggdist)
library(parallel) ## mclapply

seeds = seq(100, 9900, 100)

## Load in totalPopSpec for each seed

getTotalPopSpec = function(x) {

    gpath = paste0("/home/ben/Documents/TangledMT/Results/temporalResults/Seed_", x, "/Results/")
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

## Plot log-normal estimates over time

## Alpha
ggplot(fit, aes(g, alpha)) +
    stat_lineribbon() +
    labs(x = "Time Step", y = "Log-Series Alpha") +
    theme_classic() +
    theme(text = element_text(size = 20))

ggsave(paste0("/home/ben/Documents/TangledMT/Paper/Figures/temporalResults/alphaTime.pdf"), width = 15, height = 10)

## AlphaN
ggplot(fit, aes(g, alphaN)) +
    stat_lineribbon() +
    labs(x = "Time Step", y = "Log-Series Alpha/N") +
    theme_classic() +
    theme(text = element_text(size = 20))

ggsave(paste0("/home/ben/Documents/TangledMT/Paper/Figures/temporalResults/alphaNTime.pdf"), width = 15, height = 10)

## Average and SD
avgFit = fit %>%
    group_by(g) %>%
    summarise(avgAlpha = mean(alpha, na.rm = TRUE), 
    sdAlpha = sd(alpha, na.rm = TRUE),
    avgAlphaN = mean(alphaN, na.rm = TRUE),
    sdAlphaN = sd(alphaN, na.rm = TRUE))

## Average Alpha
ggplot(avgFit, aes(g, avgAlpha)) +
    geom_line(linewidth = 2) +
    geom_ribbon(aes(ymin = avgAlpha - sdAlpha, ymax = avgAlpha + sdAlpha), alpha = 0.25) +
    labs(x = "Time Step", y = "Log-Series Alpha") +
    theme_classic() +
    theme(text = element_text(size = 20))

ggsave(paste0("/home/ben/Documents/TangledMT/Paper/Figures/temporalResults/avgAlphaTime.pdf"), width = 15, height = 10)

## Average AlphaN
ggplot(avgFit, aes(g, avgAlphaN)) +
    geom_line(linewidth = 2) +
    geom_ribbon(aes(ymin = avgAlphaN - sdAlphaN, ymax = avgAlphaN + sdAlphaN), alpha = 0.25) +
    labs(x = "Time Step", y = "Log-Series Alpha") +
    theme_classic() +
    theme(text = element_text(size = 20))

ggsave(paste0("/home/ben/Documents/TangledMT/Paper/Figures/temporalResults/avgAlphaNTime.pdf"), width = 15, height = 10)
