#################################################
## Analyse TaNa output where it only has
## trophic interactions
#################################################

library(tidyverse)
library(ggpmisc) ## stat_poly

gpath = "/home/ben/Documents/TangledMT/Results/TNM_Output/Seed_3/Results/"
setwd(gpath)

## Load datasets

totalPop = read_delim("totalPop.txt", col_names = FALSE) %>%
    rename(g = 1, n = 2)
cellPopSpec = read_delim("cellPopSpec.txt", col_names = FALSE) %>%
    rename(g = 1, c = 2, s = 3, n = 4, M = 5, pp = 6)
totalPopSpec = read_delim("totalPopSpec.txt", col_names = FALSE) %>%
    rename(g = 1, s = 2, n = 3)
traits = read_delim("../traits.txt", col_names = FALSE) %>%
    rename(M = 1, pp = 2) %>% dplyr::select(M) %>%
    add_column(.before = "M", s = 1:nrow(.))

## Plot trait distribution of species pool
ggplot(traits, aes(log10(M), y = ..density..)) + 
    geom_histogram(fill = "grey80", col = "black") +
    theme_classic() +
    theme(text = element_text(size = 30)) +
    labs(x = "Log10(Body Mass)", y = "Density") +
    scale_y_continuous(expand = c(0, 0))

## Join trait data with totalpopspec
totalPopSpec = totalPopSpec %>% left_join(traits)

## Plot total population across the landscape over time
ggplot(totalPop, aes(g, n)) + 
    geom_line(linewidth = 2) +
    theme_classic() +
    labs(x = "Time", y = "Total Population") +
    theme(text = element_text(size = 30))

## Plot SAD including mass over species
ggplot(filter(totalPopSpec, g == max(totalPopSpec$g)) %>% slice_max(n, n = 20), aes(fct_rev(fct_reorder(as.factor(s), n)), n)) +
    geom_col() +
    theme_classic() +
    theme(text = element_text(size = 30),
    axis.text.x=element_blank()) +
    labs(x = "Species", y = "Abundance") +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 50020)) +
    geom_text(aes(label = round(log10(M), 1)), vjust = -0.25, size = 7.5)

## Does the abundance of species follow Damuth’s law? (n = M^-0.75)
ggplot(mutate(totalPopSpec, N = M^-0.75) %>% filter (g == max(totalPopSpec$g)), aes(n, N)) +
    geom_point() +
    geom_smooth(method = "lm", col = "red", linewidth = 2) +
    theme_classic() +
    labs(x = "Abundaunce in MTaNa", y = bquote("Predicted Abundaunce"~(M^-0.75))) +
    theme(text = element_text(size = 30)) +
    stat_poly_eq(use_label("eq"), size = 10) +
    stat_poly_eq(label.y = 0.9, size = 10)   

## Test Damuth’s law the other way (when logged the coefficient is the exponent)
ggplot(totalPopSpec %>% filter (g == max(totalPopSpec$g)), aes(log10(M), log10(n),)) +
    geom_point(size = 7.5, alpha = 0.5) +
    geom_smooth(method = "lm", col = "red", linewidth = 2) +
    theme_classic() +
    labs(x = "Log10(Body Mass)", y = "Log10(Abundance in MTaNa)") +
    theme(text = element_text(size = 30)) +
    stat_poly_eq(use_label("eq"), size = 10, label.x = 0.9) +
    stat_poly_eq(label.y = 0.9, label.x = 0.9, size = 10)   

## Plot abundaunce of species over time, including their mass
ggplot(filter(totalPopSpec, n > 20), aes(g, n, col = log10(M), group = log10(M))) +
    geom_line(linewidth = 2) +
    theme_classic() +
    labs(x = "Time", y = "Abundaunce", col = "Log10(Body \nMass)") +
    theme(text = element_text(size = 30)) +
    scale_colour_viridis_c()

## Biomass contribution by species over time
ggplot(filter(totalPopSpec, n > 5), aes(g, n*M, col = log10(M), group = log10(M))) +
    geom_line(linewidth = 2) +
    theme_classic() +
    labs(x = "Time", y = "Biomass", col = "Log10(Body \nMass)") +
    theme(text = element_text(size = 30)) +
    scale_colour_viridis_c()

### Change in mean mass over time, calculated as (mass*abundance)/abundance
avgMass = totalPopSpec %>% 
    mutate(totalM = n*M) %>%
    group_by(g) %>%
    summarise(avgM = sum(totalM)/sum(n))

ggplot(avgMass, aes(g, log10(avgM))) + 
    geom_line(linewidth = 1) +
    theme_classic() +
    labs(x = "Time", y = "Log10(Average Body Mass)") +
    theme(text = element_text(size = 30))

## Plot final trait distribution (not weighted by abundance)
ggplot(filter(totalPopSpec, g == max(totalPopSpec$g)), aes(log10(M), y = ..density..)) + 
    geom_histogram(fill = "grey80", col = "black") +
    theme_classic() +
    theme(text = element_text(size = 30)) +
    labs(x = "Log10(Body Mass)", y = "Density") +
    scale_y_continuous(expand = c(0, 0))

## Plot final trait distribution (weighted by abundance)
ggplot(filter(totalPopSpec, g == max(totalPopSpec$g)), aes(log10(M), y = ..density.., weight = n)) + 
    geom_histogram(fill = "grey80", col = "black") +
    theme_classic() +
    theme(text = element_text(size = 30)) +
    labs(x = "Log10(Body Mass)", y = "Density") +
    scale_y_continuous(expand = c(0, 0))

## Plot final trait distribution (weighted by biomass)
ggplot(filter(totalPopSpec, g == max(totalPopSpec$g)) %>% mutate(nM = n*M), aes(log10(M), y = ..density.., weight = nM)) + 
    geom_histogram(fill = "grey80", col = "black") +
    theme_classic() +
    theme(text = element_text(size = 30)) +
    labs(x = "Log10(Body Mass)", y = "Density") +
    scale_y_continuous(expand = c(0, 0))

## Raster figure showing which species are alive when
rasterDat = cellPopSpec %>% group_by(s) %>% summarise(N = sum(n*M)) %>% filter(N > 50) %>% distinct(s)
rasterDat = cellPopSpec %>% filter(s %in% rasterDat$s)

ggplot(mutate(rasterDat, mass = ifelse(n > 0, M, 0)), aes(g, as.factor(s), fill = log10(mass))) + 
    geom_tile() +
    theme_classic() +
    labs(x = "Time", y = "Species", fill = "Log10(Body Mass)") +
    theme(text = element_text(size = 30),
    axis.text.y = element_blank()) + 
    scale_fill_viridis_c()

### Analyse consumption data
consumption = read_delim(paste0(gpath, "consumptionRate.txt"), col_names = FALSE) %>%
    rename(g = 1, c = 2, Si = 3, Mi = 4, Ni = 5, Sj = 6, Mj = 7, Nj = 8, aij = 9, hij = 10, Jij = 11) %>%
    mutate(eNjJij = 0.5*Nj*Jij, NiJij = Ni*Jij)

## NjJij = consumption rate of i feeding on j when i is the focal individual/species (population of i is one)
## NiJij = consumption rate of i feeding on j, when j is the focal individual/species (population of j is one)
## Confirmed in MTaNa that this is correct

## Why are some species so abundant?
cellPopSpec %>% filter(g == max(totalPopSpec$g) & c == 1) %>% arrange(-n)

consumption %>% filter(g == max(consumption$g) & Si == 229 & c == 1) %>% arrange(-eNjJij)

consumption %>% filter(g == max(consumption$g) & Sj == 229 & c == 1) %>% arrange(-NiJij) %>% print(w = 100)

###############################
## Food Web
###############################

library(igraph)
library(ggnetwork)
library(sna)

n = network(rgraph(10, tprob = 0.2), directed = FALSE)

web = read_delim("consumptionRate.txt", col_names = FALSE) %>%
    rename(g = 1, c = 2, Si = 3, Mi = 4, Ni = 5, Sj = 6, Mj = 7, Nj = 8, C = 9)

## Change to be names of species i and j in first two columns
web = web %>% relocate(Si, Sj, .before = g)

test = web %>% filter(g == 100000)
testMat = test %>% dplyr::select(Si, Sj, C) %>% pivot_wider(names_from = Sj, values_from = C, values_fill = 0) %>%
    column_to_rownames("Si")

ggplot(ggnetwork(test), aes(x, y, xend = xend, yend = yend)) + 
    geom_edges(aes(linewidth = log10(C/Mi)), color = "grey50", curvature = 0.1,
    arrow = arrow(length = unit(6, "pt"), type = "closed")) +
    geom_nodes(aes(size = log10(Nj*Mj))) +
    geom_nodetext(aes(label = round(log10(Mi), 2)), fontface = "bold", col = "white") +
    geom_edgetext(aes(label = round(C/Mi, 2)), color = "white", fill = "grey25",
    box.padding = unit(1, "lines")) +
    theme_blank()
