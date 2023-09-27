#################################################
## Analyse best distribution to fit
## SAD at each time step
## for each temporal model
## output : best distribution for each time step
## and calculate % of time steps for which
## each distribution is the best
#################################################

library(tidyverse)
library(fitdistrplus) ## fit distributions to data
library(sads) ## fit logseries dist


seeds = seq(100, 10000, 100)

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

testDists = function(x) {

    seed = x$seed[[1]]
    g = x$g[[1]]

    lnormFit = tryCatch(fitdist(x$n, "lnorm"), error = function(e) {NA})
    lsFit = tryCatch(fitsad(x$n, "ls"), error = function(e) {NA})
    unFit = tryCatch(fitdist(x$n, "unif"), error = function(e) {NA})
    nFit = tryCatch(fitdist(x$n, "norm"), error = function(e) {NA})

    out = data.frame(lnorm = tryCatch(logLik(lnormFit), error = function(e) {NA}),
        ls = tryCatch(logLik(lsFit), error = function(e) {NA}),
        unif = tryCatch(logLik(unFit), error = function(e) {NA}),
        norm = tryCatch(logLik(nFit), error = function(e) {NA})) %>%
        mutate(.before = 1, seed, g)

    return(out)
}

dists = totalPopSpec %>%
    group_by(seed, g) %>%
    group_map(~testDists(.), .keep = TRUE) %>%
    bind_rows()

write_csv(dists, paste0(gpath, "../../../../Results/temporalResults/SADDistribution.csv"))

## Calculate best distribution for each time step in each seed
bestDist = dists %>%
    pivot_longer(cols = c(-seed, -g), names_to = "dist", values_to = "logLik") %>%
    group_by(seed, g) %>%
    slice_max(logLik) %>%
    mutate(dist = case_when(dist == "lnorm" ~ "Log-Normal", dist == "ls" ~ "Log-Series", dist == "unif" ~ "Uniform", dist == "norm" ~ "Normal"))

plotDist = bestDist %>%
    group_by(dist) %>%
    summarise(n = n()) %>%
    mutate(p = n/sum(n))

## Plot the bestDist
ggplot(plotDist, aes(reorder(dist, -p), p)) +
    geom_col(col = "black", size = 1.5, fill = "grey70") +
    theme_classic() +
    labs(x = "Best Distribution (Log-Likelihood)", y = "Proportion") +
    theme(text = element_text(size = 30)) +
    scale_y_continuous(expand = c(0, 0), limits = c(0,1.1)) +
    geom_text(aes(dist, y = p, label = round(p, 2)), vjust = -0.25, size = 8)

ggsave(paste0(gpath, "../../../../Paper/Figures/temporalResults/bestSADDistribution.pdf"), width = 15, height = 10)


