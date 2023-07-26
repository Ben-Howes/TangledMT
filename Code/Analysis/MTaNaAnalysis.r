#################################################
## Analyse TaNa output where it only has
## trophic interactions
#################################################

library(tidyverse)
library(ggpmisc) ## stat_poly
library(lemon) ##facet_rep_wrap

gpath = "/home/ben/Documents/TangledMT/Results/TNM_Output/Seed_2/Results/"
setwd(gpath)

## Load datasets

totalPop = read_delim("totalPop.txt", col_names = FALSE) %>%
    rename(g = 1, n = 2)
cellPopSpec = read_delim("cellPopSpec.txt", col_names = FALSE) %>%
    rename(g = 1, c = 2, s = 3, n = 4, M = 5, pp = 6)
totalPopSpec = read_delim("totalPopSpec.txt", col_names = FALSE) %>%
    rename(g = 1, s = 2, n = 3)
traits = read_delim("../traits.txt", col_names = FALSE) %>%
    rename(M = 1, pp = 2) %>%
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

# ggsave(paste0(gpath, "../../../../Paper/Figures/InitialMTaNaFauna/totalPop.pdf"), width = 18, height = 10)

## Plot SAD but with mass
ggplot(filter(totalPopSpec, g == max(totalPopSpec$g)) %>% slice_max(n, n = 100), aes(log10(M), log10(n), fill = as.factor(pp))) +
    geom_col(width = 0.01) +
    theme_classic() +
    theme(text = element_text(size = 30)) +
    labs(x = "Log10(Mass)", y = "Log10(Abundance", fill = "Primary\nProducer") +
    scale_y_continuous(expand = c(0, 0)) +
    scale_fill_viridis_d()

# ggsave(paste0(gpath, "../../../../Paper/Figures/InitialMTaNaFauna/SAD.pdf"), width = 18, height = 10)

## Test Damuth’s law the other way (when logged the coefficient is the exponent)
ggplot(totalPopSpec %>% filter (g == max(totalPopSpec$g) & n > 1), aes(log10(M), log10(n), col = as.factor(pp))) +
    geom_point(size = 7.5, alpha = 0.5) +
    geom_smooth(method = "lm", linewidth = 2) +
    theme_classic() +
    labs(x = "Log10(Body Mass)", y = "Log10(Abundance in MTaNa)",
    col = "Primary\nProducer") +
    theme(text = element_text(size = 30)) +
    stat_poly_eq(use_label("eq"), size = 10, label.x = 0.9) +
    scale_colour_viridis_d(end = 0.7)

checkDamuth = function(x) {
    x = x %>% filter(pp == 0)
    if(nrow(x) > 4) {
        mod = lm(log10(n) ~ log10(M), data = x)
        slope = coef(mod)[[2]]
        out = x %>% mutate(dam = slope)
        return(out)
        }
}

damuth = cellPopSpec %>% 
    filter(g %% 100000 == 0) %>%
    group_split(g) %>%
    lapply(., checkDamuth) %>% 
    bind_rows()

ggplot(damuth, aes(g, dam)) + geom_line()

## Plot abundance of species over time, including their mass
ggplot(filter(totalPopSpec), aes(g, log10(n), col = log10(M), group = interaction(log10(M), as.factor(pp)), linetype = as.factor(pp))) +
    geom_line(linewidth = 2) +
    theme_classic() +
    labs(x = "Time", y = "Log10(Abundaunce)", col = "Log10(Body \nMass)", 
    linetype = "Primary Producer") +
    theme(text = element_text(size = 30)) +
    scale_colour_viridis_c()

# ggsave(paste0(gpath, "../../../../Paper/Figures/InitialMTaNaFauna/SpeciesAbundanceTime.pdf"), width = 18, height = 10)

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

## Raster figure showing which species are alive when
rasterDat = cellPopSpec %>% group_by(s) %>% summarise(N = sum(n*M)) %>% filter(N > 50) %>% distinct(s)
rasterDat = cellPopSpec %>% filter(s %in% rasterDat$s)

ggplot(mutate(rasterDat, mass = ifelse(n > 0, M, 0)), aes(g, as.factor(s), fill = log10(mass))) + 
    geom_tile() +
    facet_rep_wrap(. ~ pp) +
    theme_classic() +
    labs(x = "Time", y = "Species", fill = "Log10(Body Mass)") +
    theme(text = element_text(size = 30),
    axis.text.y = element_blank()) + 
    scale_fill_viridis_c()

### Analyse consumption data
consumption = read_delim(paste0(gpath, "consumptionRate.txt"), col_names = FALSE) %>%
    rename(g = 1, c = 2, Si = 3, Mi = 4, Ni = 5, Sj = 6, Mj = 7, Nj = 8, aij = 9, Aij = 10, hij = 11, Jij = 12) %>%
    mutate(NjJij = Nj*Jij, eNjJij = 0.5*Nj*Jij, NiJij = Ni*Jij)

## NjJij = consumption rate of i feeding on j when i is the focal individual/species (population of i is one)
## NiJij = consumption rate of i feeding on j, when j is the focal individual/species (population of j is one)
## Confirmed in MTaNa that this is correct

## Calculate H for each time step for each species
gain = consumption %>% group_by(g, c, Si) %>% summarise(gain = sum(eNjJij))
loss = consumption %>% group_by(g, c, Sj) %>% summarise(loss = sum(NiJij))
z = consumption %>% distinct(Si, Mi) %>% mutate(z = (4.15*(10^-8))*Mi^0.75)

calculateSearchRate = function(mi, mj, T) {

    V0 = 0.33
    D0 = 1.62
    P0 = 1 ## Temperature constant (vary by taxa)
    k = 8.6173*(10^-5) ## Boltzmann constant
    T0 = 293.15 ## 0 celsius in Kelvin
    E = 0 ## Activation energy

    aij = 2*(V0)*(D0)*(mi^(0.63))*(exp(1)^(-E/(k*(T + T0))))

    return(aij)

}

density = consumption %>% group_by(g, c, Si) %>% distinct(Si, Mi, Ni) %>% mutate(NiJii = 0.1*Mi*Ni*calculateSearchRate(Mi, Mi, 0))

joined = left_join(gain, loss, by = join_by(Si == Sj, g, c)) %>% left_join(z) %>% left_join(density) %>%
    mutate(H = gain - loss - NiJii - z) %>% 
    mutate(HM = H/Mi)

joined = joined %>% left_join(dplyr::select(traits, s, pp), by = join_by(Si == s)) %>%
    mutate(gain = ifelse(pp == 1, Mi*(1 - (Ni*Mi)/(10*(Mi^0.25))), gain),
    H = ifelse(pp == 1, gain - loss - z, H), HM = (H/Mi), pOff = (1 / (1 + exp(-(1/(Mi^0.25))*(HM - 0.5)))))

ggplot(filter(joined, Si == 219), aes(g, Ni, col = cut(pOff, c(-Inf, 0.15, Inf)))) +
    geom_point(size = 2) +
    scale_colour_manual(name = "pOff", values = c("(-Inf,0.15]" = "red",
                                  "(0.15, Inf]" = "blue"),
                                  labels = c("<= 0.15", "0.15 <")) +
    theme_classic() +
    theme(text = element_text(size = 30))

###############################
## Food Web
###############################

library(diagram)

web = test1

## Change to be names of species i and j in first two columns
web = web %>% relocate(Si, Sj, .before = g)

testMat = web %>% mutate(eNjJij = eNjJij/Mi) %>% dplyr::select(Si, Sj, eNjJij) %>% pivot_wider(names_from = Sj, values_from = eNjJij, values_fill = 0) %>%
    column_to_rownames("Si")

testMat = testMat[, order(as.numeric(colnames(testMat)))]
testMat = testMat[order(as.numeric(row.names(testMat))), ]
testMat[testMat < 1] = 0 

plotmat(testMat, relsize = 0.8)
