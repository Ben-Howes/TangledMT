#################################################
## Calculate total population and
## richness over time
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

totalPop = totalPopSpec %>% 
    group_by(seed, g) %>%
    summarise(n = sum(n))

totalRich = totalPopSpec %>% 
    group_by(seed, g) %>%
    summarise(n = n())

## Plot over time

## Total Population
ggplot(totalPop, aes(g, n)) +
    stat_lineribbon() +
    labs(x = "Time Step", y = "Total Population") +
    theme_classic() +
    theme(text = element_text(size = 20))

ggsave(paste0("/home/ben/Documents/TangledMT/Paper/Figures/temporalResults/totalPopTime.pdf"), width = 15, height = 10)

## Total Richness
ggplot(totalRich, aes(g, n)) +
    stat_lineribbon() +
    labs(x = "Time Step", y = "Total Richness") +
    theme_classic() +
    theme(text = element_text(size = 20))

ggsave(paste0("/home/ben/Documents/TangledMT/Paper/Figures/temporalResults/totalRichTime.pdf"), width = 15, height = 10)
