#################################################
## Combine damuth law over time
## across all seeds and each time step
## to see if it changes with climate change
#################################################

library(tidyverse)
library(ggdist)
library(parallel) ## mclapply

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

 ## Test damuth law over time
checkDamuth = function(x) {
    x = x %>% filter(pp == 0 & n > 4)
    if(nrow(x) > 4) {
        mod = lm(log10(n) ~ log10(M), data = x)
        slope = coef(mod)[[2]]
        out = x %>% distinct(seed, g) %>% mutate(dam = slope)
        return(out)
    }
}

damuth = totalPopSpec %>% 
    group_split(seed, g) %>%
    mclapply(., checkDamuth, mc.cores = 8) %>% 
    bind_rows()

ggplot() +
    geom_rect(aes(xmin = 0, xmax = tClimate, ymin = -Inf, ymax = Inf), fill = "lightblue", alpha = 0.5) +
    geom_rect(aes(xmin = tClimate, xmax = tMax, ymin = -Inf, ymax = Inf), fill = "red", alpha = 0.25) +
    stat_lineribbon(data = damuth, aes(g, dam)) +
    geom_vline(xintercept = tClimate, linetype = "dashed", linewidth = 2) +
    geom_hline(yintercept = -0.75, linetype = "dashed", linewidth = 2, colour = "red") +
    labs(x = "Time Step", y = "Mass-Abundance Exponent") +
    theme_classic() +
    theme(text = element_text(size = 20)) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_continuous(expand = c(0, 0))

ggsave(paste0("/home/ben/Documents/TangledMT/Paper/Figures/climateChange/climateDamuth.pdf"), width = 15, height = 10)
ggsave(paste0("/home/ben/Documents/TangledMT/Paper/Figures/climateChange/climateDamuth.png"), width = 15, height = 10)

avgDam = damuth %>%
    group_by(g) %>%
    summarise(avgDam = mean(dam, na.rm = TRUE), 
    sd = sd(dam, na.rm = TRUE))

ggplot() +
    geom_rect(aes(xmin = 0, xmax = tClimate, ymin = -Inf, ymax = Inf), fill = "lightblue", alpha = 0.5) +
    geom_rect(aes(xmin = tClimate, xmax = tMax, ymin = -Inf, ymax = Inf), fill = "red", alpha = 0.25) +
    geom_line(data = avgDam, aes(g, avgDam), linewidth = 2) +
    geom_ribbon(data = avgDam, aes(x = g, ymin = avgDam - sd, ymax = avgDam + sd), alpha = 0.5) +
    geom_vline(xintercept = tClimate, linetype = "dashed", linewidth = 2) +
    geom_hline(yintercept = -0.75, linetype = "dashed", linewidth = 2, colour = "red") +
    labs(x = "Time Step", y = "Mass-Abundance Exponent") +
    theme_classic() +
    theme(text = element_text(size = 20)) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_continuous(expand = c(0, 0))

ggsave(paste0("/home/ben/Documents/TangledMT/Paper/Figures/climateChange/climateAvgDamuth.pdf"), width = 15, height = 10)
ggsave(paste0("/home/ben/Documents/TangledMT/Paper/Figures/climateChange/climateAvgDamuth.png"), width = 15, height = 10)

## Try some stats
## Test if damuthis significantly different before and after climate change
## We need to remove the "burn in" values, and I need to figure out exactly when this happens
## I think going from 2,500,000 is a safe conservative place to start

modData = damuth %>% 
    filter(g >= 25000000) %>%
    mutate(climate = ifelse(g < tClimate, 0, 1))

mod = lm(dam ~ climate, data = modData)

ggplot(modData, aes(as.factor(climate), dam)) +
    geom_boxplot(aes(fill = as.factor(climate)), alpha = c(0.5, 0.25)) +
    theme_classic() +
    labs(x = "Climate Change", y = "Mass-Abundance Exponent", fill = NULL) +
    theme(legend.position = "none", text = element_text(size = 20)) +
    scale_x_discrete(labels = c("Before", "After")) +
    scale_fill_manual(values = c("lightblue", "red"))

ggsave(paste0("/home/ben/Documents/TangledMT/Paper/Figures/climateChange/damuthDiff.pdf"), width = 15, height = 10)
ggsave(paste0("/home/ben/Documents/TangledMT/Paper/Figures/climateChange/damuthDiff.png"), width = 15, height = 10)
