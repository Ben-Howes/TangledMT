############################################
## Script to analyse outputs of MTaNa
## using different constants
############################################

library(tidyverse)
library(ggpmisc) ## stat_poly
library(lemon) ##facet_rep_wrap
library(parallel) ## mclapply

gpath = "/home/ben/Documents/TangledMT/Results/MTaNaConstants/"
setwd(gpath)

## We want to analyse simple metrics between models run with different constants
## I plan to test differences in:
## 1) Total Abundance, 2) Total richness
## 3) Body mass distribution 4) Dalmuth rule constant (per gen)

## Get paths to files from each run
paths = list.files(gpath)

## Read in parameters, seed, and cellPopInd for each run - join cellPopInd with parameters

cellPopSpec = function(path) {
    params = read_table(paste0(gpath, path, "/Parameters.txt"), col_names = F) %>%
        dplyr::select(X2, X3) %>% rename("constant" = 1, "value" = 2) %>%
        pivot_wider(names_from = "constant", values_from = "value")
    cellPopSpec = read_table(paste0(gpath, path, "/Results/cellPopSpec.txt"), col_names = F) %>%
        rename("g" = 1, "c" = 2, "s" = 3, "n" = 4, "m" = 5, "pp" = 6) %>%
        mutate(.before = "g", seed = path)
    cellPopSpec = cbind(cellPopSpec, params)

    return(cellPopSpec)
}

cellPopSpec = mclapply(paths, cellPopSpec, mc.cores = 6) %>% 
    bind_rows() %>%
    filter(g > 500000)

##################################
## Total Population
##################################

totalPop = cellPopSpec %>% group_by(seed, g, r0, K0, I0) %>% 
    filter(pp == 0) %>% summarise(totalPop = sum(n))

ggplot(totalPop, aes(g, log10(totalPop), col = as.factor(I0), group = as.factor(seed))) +
    facet_rep_grid(as.factor(r0) ~ as.factor(K0)) +
    geom_line(size = 1) +
    scale_colour_viridis_d() + 
    labs(x = "Time", y = "Log10(Abundance)", col = "Intraspecific\nCompetition (I0)",
    title = "Rows = Primary Producer Growth Rate (r0)\nCols = Primary Producer Carrying Capacity (K0)") +
    scale_x_continuous(n.breaks = 3) + 
    theme_classic() + 
    theme(text = element_text(size = 30))

##################################
## Total Species Richness
##################################

totalRich = cellPopSpec %>% group_by(seed, g, r0, K0, I0) %>% 
    filter(pp == 0) %>%
    summarise(totalRich = n_distinct(s))

ggplot(totalRich, aes(g, log10(totalRich), col = as.factor(I0), group = as.factor(seed))) +
    facet_rep_grid(as.factor(r0) ~ as.factor(K0)) +
    geom_line(size = 1) +
    scale_colour_viridis_d() + 
    labs(x = "Time", y = "Log10(Richness)", col = "Intraspecific\nCompetition (I0)",
    title = "nRows = Primary Producer Growth Rate (r0)\nCols = Primary Producer Carrying Capacity (K0)") +
    scale_x_continuous(n.breaks = 3) + 
    theme_classic() + 
    theme(text = element_text(size = 30))

##################################
## Average body mass
##################################

avgMass = cellPopSpec %>% 
    mutate(biomass = m*n) %>% 
    group_by(seed, g, r0, K0, I0) %>% 
    filter(pp == 0) %>%
    summarise(avgMass = sum(biomass)/n())

ggplot(avgMass, aes(g, log10(avgMass), col = as.factor(I0), group = as.factor(seed))) +
    facet_rep_grid(as.factor(r0) ~ as.factor(K0)) +
    geom_line(size = 1) +
    scale_colour_viridis_d() + 
    labs(x = "Time", y = "Log10(Average Mass)", col = "Intraspecific\nCompetition (I0)",
    title = "Rows = Primary Producer Growth Rate (r0)\nCols = Primary Producer Carrying Capacity (K0)") +
    scale_x_continuous(n.breaks = 3) + 
    theme_classic() + 
    theme(text = element_text(size = 30))

##################################
## Does it follow Damuth's law?
##################################

## Run an lm for each time point (every 1000 time steps)
## of log10(n) as explained by log10(m)
## save the slope of log10(m) which should be -0.75
## Only do this for non-producers since producers follow damuth's law
## based on the equations we've set

checkDamuth = function(x) {
    x = x %>% filter(pp == 0)
    if(nrow(x) > 4) {
        mod = lm(log10(n) ~ log10(m), data = x)
        slope = coef(mod)[[2]]
        out = x %>% mutate(dam = slope)
        return(out)
        }
}

damuth = cellPopSpec %>% 
    filter(g %% 100000 == 0) %>%
    group_split(seed, g, r0, K0, I0) %>%
    lapply(., checkDamuth) %>% 
    bind_rows()

ggplot(damuth, aes(g, dam, col = as.factor(I0), group = as.factor(seed))) +
    facet_rep_grid(as.factor(r0) ~ as.factor(K0)) +
    geom_line(linewidth = 1) +
    scale_colour_viridis_d() + 
    labs(x = "Time", y = "Slope of Log10(Abundance) ~ Log10(Mass)", col = "Intraspecific\nCompetition (I0)",
    title = "Rows = Primary Producer Growth Rate (r0)\nCols = Primary Producer Carrying Capacity (K0)") +
    scale_x_continuous(n.breaks = 3) + 
    theme_classic() + 
    theme(text = element_text(size = 30)) +
    geom_hline(yintercept = -0.75, linetype = "dashed")
