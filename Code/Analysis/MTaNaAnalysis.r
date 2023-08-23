#################################################
## Analyse TaNa output where it only has
## trophic interactions
#################################################

library(tidyverse)
library(ggpmisc) ## stat_poly
library(lemon) ##facet_rep_wrap
library(colourvalues) ## colour igraph vertices
library(fields) ## add legend to igraph

seeds = seq(100, 500, 100)

for (x in seeds) {

    gpath = paste0("/home/ben/Documents/TangledMT/Results/testParallel/Seed_", x, "/Results/")
    setwd(gpath)

    dir.create(paste0(gpath, "../../../../Paper/Figures/testParallel/Seed_", x), showWarnings = FALSE)

    ## Load datasets
    totalPop = read_delim("totalPop.txt", col_names = FALSE) %>%
        rename(g = 1, n = 2)
    totalRich = read_delim("totalRich.txt", col_names = FALSE) %>%
        rename(g = 1, n = 2)
    cellPopSpec = read_delim("cellPopSpec.txt", col_names = FALSE) %>%
        rename(g = 1, c = 2, s = 3, n = 4, biomass = 5, M = 6, pp = 7)
    totalPopSpec = read_delim("totalPopSpec.txt", col_names = FALSE) %>%
        rename(g = 1, s = 2, n = 3)
    traits = read_delim("../traits.txt", col_names = FALSE) %>%
        rename(M = 1, pp = 2) %>%
        add_column(.before = "M", s = 1:nrow(.))

    ## Plot trait distribution of species pool
    # ggplot(traits, aes(log10(M), y = ..density..)) + 
    #     geom_histogram(fill = "grey80", col = "black") +
    #     theme_classic() +
    #     theme(text = element_text(size = 30)) +
    #     labs(x = "Log10(Body Mass)", y = "Density") +
    #     scale_y_continuous(expand = c(0, 0))

    ## Join trait data with totalpopspec
    totalPopSpec = totalPopSpec %>% left_join(traits)

    ## Plot total population across the landscape over time
    ggplot(totalPop, aes(g, n)) + 
        geom_line(linewidth = 2) +
        theme_classic() +
        labs(x = "Time", y = "Total Population") +
        theme(text = element_text(size = 30))

    ggsave(paste0(gpath, "../../../../Paper/Figures/testParallel/Seed_", x, "/totalPop_", x, ".pdf"), width = 15, height = 10)

    ## Plot total richness across the landscape over time
    ggplot(totalRich, aes(g, n)) + 
        geom_line(linewidth = 2) +
        theme_classic() +
        labs(x = "Time", y = "Species Richness") +
        theme(text = element_text(size = 30))

    ggsave(paste0(gpath, "../../../../Paper/Figures/testParallel/Seed_", x, "/totalRich_", x, ".pdf"), width = 15, height = 10)

    ## Plot SAD but with mass
    ggplot(filter(totalPopSpec, g == max(totalPopSpec$g)), aes(log10(M), log10(n), fill = as.factor(pp))) +
        geom_col(width = 0.01) +
        theme_classic() +
        theme(text = element_text(size = 30)) +
        labs(x = "Log10(Mass)", y = "Log10(Abundance", fill = NULL) +
        scale_y_continuous(expand = c(0, 0)) +
        scale_fill_viridis_d(labels = c("Consumer", "Primary Producer"))

    ggsave(paste0(gpath, "../../../../Paper/Figures/testParallel/Seed_", x, "/SAD_", x, ".pdf"), width = 15, height = 10)

    ## Test Damuth’s law the other way (when logged the coefficient is the exponent)
    ggplot(totalPopSpec %>% filter (g == max(totalPopSpec$g) & n > 4), aes(log10(M), log10(n), col = as.factor(pp))) +
        geom_point(size = 7.5, alpha = 0.5) +
        geom_smooth(method = "lm", linewidth = 2) +
        theme_classic() +
        labs(x = "Log10(Body Mass)", y = "Log10(Abundance in MTaNa)", 
        col = NULL) +
        theme(text = element_text(size = 30)) +
        stat_poly_eq(use_label("eq"), size = 10, label.x = 0.9) +
        scale_colour_viridis_d(end = 0.7, labels = c("Consumer", "Primary Producer"))

    ggsave(paste0(gpath, "../../../../Paper/Figures/testParallel/Seed_", x, "/damuthLinear_", x, ".pdf"), width = 15, height = 10)

    checkDamuth = function(x) {
        x = x %>% filter(pp == 0 & n > 4)
        if(nrow(x) > 4) {
            mod = lm(log10(n) ~ log10(M), data = x)
            slope = coef(mod)[[2]]
            out = x %>% distinct(g, c) %>% mutate(dam = slope)
            return(out)
        }
    }

    damuth = cellPopSpec %>% 
        filter(g %% 100000 == 0) %>%
        group_split(g) %>%
        lapply(., checkDamuth) %>% 
        bind_rows()

    ggplot(damuth, aes(g, dam)) + 
        geom_line(linewidth = 1) +
        theme_classic() +
        geom_hline(yintercept = -0.75, linetype = "dashed", linewidth = 1) +
        theme(text = element_text(size = 30)) +
        labs(x = "Time", y = "Damuth's Exponent")

    ggsave(paste0(gpath, "../../../../Paper/Figures/testParallel/Seed_", x, "/damuthTime_", x, ".pdf"), width = 15, height = 10)

    ## Plot abundance of species over time, including their mass
    ggplot(filter(totalPopSpec), aes(g, log10(n), col = log10(M), group = interaction(log10(M), as.factor(pp)), linetype = as.factor(pp))) +
        geom_line(linewidth = 2) +
        theme_classic() +
        labs(x = "Time", y = "Log10(Abundaunce)", col = "Log10(Body \nMass)", 
        linetype = "Primary Producer") +
        theme(text = element_text(size = 30)) +
        scale_colour_viridis_c()

    # ggsave(paste0(gpath, "../../../../Paper/Figures/InitialMTaNaFauna/SpeciesAbundanceTime.pdf"), width = 18, height = 10)

    ################################
    ## Reproduction Data
    ################################

    # reproduction = read_delim(paste0(gpath, "reproduction.txt"), col_names = FALSE) %>% 
    #     rename(g = 1, c = 2, s = 3, m = 4, H = 5, pOffMax = 6, pOff = 7, pDeath = 8) %>%
    #     mutate(HM = H/m, .after = H) %>%
    #     left_join(traits)

    # ## Distribution of H across all time steps for all species
    # ggplot(reproduction, aes(HM)) + 
    #     geom_density(linewidth = 2) +
    #     theme_classic() +
    #     labs(x = "H'", y = "Density") +
    #     theme(text = element_text(size = 30))

    # ## Proportion of maxPoff achieved per time step for each species
    # ggplot(reproduction, aes(g, pOff/pOffMax, col = log10(m), group = log10(m))) +
    #     geom_line(linewidth = 1) +
    #     theme_classic() +
    #     labs(x = "Time", y = "Max Probablity of Reproduction %") +
    #     theme(text = element_text(size = 30)) +
    #     scale_colour_viridis_c()

    # ## pOff - Pdeath for each species in each timestep
    # ggplot(reproduction, aes(g, pOff - pDeath, col = log10(m), group = log10(m))) +
    #     geom_line(linewidth = 1) +
    #     theme_classic() +
    #     labs(x = "Time", y = "(Probability of Reproduction) - (Probability of Death)") +
    #     theme(text = element_text(size = 30)) +
    #     scale_colour_viridis_c() +
    #     geom_hline(yintercept = 0, linetype = "dashed", linewidth = 2)

    # ## Probability of reproduction for each species in each time step 
    # ## with line for maxPoff and pDeath
    # ggplot(reproduction, aes(log10(m), pOff)) + 
    #     geom_point(size = 2.5, shape = 1) + 
    #     geom_line(aes(log10(m), pOffMax), linewidth = 2, col = "blue") +
    #     geom_line(aes(log10(m), pDeath), linewidth = 2, col = "red") +
    #     theme_classic() +
    #     labs(x = "Log10(Body Mass)", y = "Probability of Reproduction",
    #     title = "Probability of Reproduction\nBlue Line is Maximum Probability of Reproduction\nRed Line is Probability of Death") +
    #     theme(text = element_text(size = 30))

    # ## Maximum average number of offspring produced by an individual will be pR/pZ, wheree pR is it it's max value of 1/G
    # reproduction = reproduction %>% mutate(avgOffspring = pOff/pDeath, maxOffspring = pOffMax/pDeath)

    # ggplot(reproduction, aes(log10(m), maxOffspring)) +
    #     geom_point(size = 2.5, shape = 1) + 
    #     theme_classic() +
    #     labs(x = "Log10(Body Mass)", y = "Maximum Number of Offspring in a lifetime") +
    #     theme(text = element_text(size = 30))

    # ggplot(reproduction, aes(log10(m), avgOffspring)) +
    #     geom_point(size = 2.5, shape = 1) + 
    #     theme_classic() +
    #     labs(x = "Log10(Body Mass)", y = "Average Number of Offspring in a lifetime") +
    #     theme(text = element_text(size = 30))

    ################################
    ## Consumption Rate
    ################################

    ### Analyse consumption data
    consumption = read_delim(paste0(gpath, "consumptionRate.txt"), col_names = FALSE) %>%
        rename(g = 1, c = 2, Si = 3, Mi = 4, Ni = 5, Sj = 6, Mj = 7, Nj = 8, Jij = 9) %>%
        mutate(Jij = as.numeric(ifelse(Jij == "0 0 0 0", 0, Jij))) %>%
        mutate(NjJij = Nj*Jij, eNjJij = 0.5*Nj*Jij, NiJij = Ni*Jij)

    ## NjJij = consumption rate of i feeding on j when i is the focal individual/species (population of i is one)
    ## NiJij = consumption rate of i feeding on j, when j is the focal individual/species (population of j is one)
    ## Confirmed in MTaNa that this is correct

    ## Now we want to make a food web from the consumption rate
    ## Since all individuals can feed on each other, we need to use some sort of threshold
    ## which is based on the mass of the consumer - probably a percentage of the consumer mass

    library(igraph)
    library(fluxweb)

    ## Load in example data
    ## and function for plotfw
    # load(paste0(gpath, "../../../../../../Downloads/BalticFoodWeb/BalticFW.Rdata")) 
    source(paste0(gpath, "../../../../Code/Analysis/Rscripts/networkPlotScripts.r"))

    ## Filter dataset to a specific time step, and summarise
    ## over cells to get a total consumption rate of each Si on each Sj
    test = consumption %>% 
        filter(g == 50000000) %>%
        group_by(Si, Ni, Mi, Sj, Mj, Nj) %>%
        summarise(Ni = sum(Ni), Nj = sum(Nj), eNjJij = sum(eNjJij)) %>%
        ungroup()

    ## Remove 0 interactions (basal species), and those below a threshold (using non-zero species mean as the current threshold)
    ## and relocate so Si and Sj are first as needed in graph_from_data_frame function
    cDF = test %>%
        mutate(gain = eNjJij/Mi) %>%
        filter(gain > 0) %>%
        filter(gain > mean(gain)) %>%
        dplyr::select(Sj, Si, Nj, Ni, gain) %>%
        relocate(Sj, Si, Nj, Ni, gain) %>%
        mutate(gainPalette = alpha("black", gain/max(gain)))

    ## Make another data frame with species names, mass and abundances
    mDF = totalPopSpec %>% 
        filter(g == 50000000) %>%
        dplyr::select(s, M, n, pp) %>%
        rename(m = 2) %>%
        mutate(logM = log10(m)) %>%
        mutate(massPalette = colour_values(logM, palette = "viridis"))

    ## Make igraph object from data frames
    ## using vertices = mDF means we get the metadata into the object as well
    ## so can easily change colour/size based on abundance or mass
    testG = graph_from_data_frame(cDF, directed = TRUE, vertices = mDF)

    ## Scale biomass size
    minN = 3 ## minimum size of vertex
    scaleN = 10 ## maximum size of vertex
    maxN = max(V(testG)$n)
    sizeN = (V(testG)$n/maxN)*scaleN + minN

    ## Change colour based on log10(mass)
    col = arrange(mDF, logM)

    pdf(file = paste0(gpath, "../../../../Paper/Figures/testParallel/Seed_", x, "/foodWeb_", x, ".pdf"), width = 16, height = 16) 
    plotfw(testG, edge.arrow.size = 0.3, size = sizeN, vertex.color = V(testG)$massPalette,
        edge.color = E(testG)$gainPalette, edge.width = 4)
    image.plot(legend.only=T, zlim=range(mDF$logM), col = col$massPalette)
    dev.off()

    ###############
    ## HEATMAP
    ###############

    # ## Plot heatmap of interactions
    # mat = test %>% 
    #     mutate(gain = log10(eNjJij/Mi)) %>%
    #     dplyr::select(Si, Sj, gain) %>%
    #     pivot_wider(names_from = Sj, values_from = gain) %>%
    #     column_to_rownames("Si") %>%
    #     as.matrix()

    # heatmap(mat, Rowv = NA, Colv = NA, scale="none")
}
