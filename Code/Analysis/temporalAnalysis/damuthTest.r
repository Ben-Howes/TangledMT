#################################################
## Combine damuth law over time
## across all seeds and each time step
## to get a general idea of whether
## communities follow mass distributions
#################################################

library(tidyverse)
library(ggdist)
library(parallel) ## mclapply

seeds = seq(100, 9900, 100)

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

ggplot(damuth, aes(g, dam)) +
    stat_lineribbon() +
    geom_hline(yintercept = -0.75, linetype = "dashed", linewidth = 2, colour = "red") +
    labs(x = "Time Step", y = "Mass-Abundance Exponent") +
    theme_classic() +
    theme(text = element_text(size = 20)) +
    ylim(-1, 1)

ggsave(paste0("/home/ben/Documents/TangledMT/Paper/Figures/temporalResults/overallDamuth.pdf"), width = 15, height = 10)

avgDam = damuth %>%
    group_by(g) %>%
    summarise(avgDam = mean(dam, na.rm = TRUE), 
    sd = sd(dam, na.rm = TRUE))

ggplot(avgDam, aes(g, avgDam)) +
    geom_line(linewidth = 2) +
    geom_ribbon(aes(ymin = avgDam - sd, ymax = avgDam + sd), alpha = 0.25) +
    geom_hline(yintercept = -0.75, linetype = "dashed", linewidth = 2, colour = "red") +
    labs(x = "Time Step", y = "Mass-Abundance Exponent") +
    theme_classic() +
    theme(text = element_text(size = 20)) +
    ylim(-1,1)

ggsave(paste0("/home/ben/Documents/TangledMT/Paper/Figures/temporalResults/avgDamuth.pdf"), width = 15, height = 10)
