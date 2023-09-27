#################################################
## Analyse TaNa output where it only has
## trophic interactions
#################################################

library(tidyverse)
library(ggpmisc) ## stat_poly
library(lemon) ##facet_rep_wrap
library(colourvalues) ## colour igraph vertices
library(fields) ## add legend to igraph

seeds = seq(100, 10000, 100)

for (x in seeds) {

    gpath = paste0("/home/ben/Documents/TangledMT/Results/climateChange/Seed_", x, "/Results/")
    setwd(gpath)

    dir.create(paste0(gpath, "../../../../Paper/Figures/climateChange/Seed_", x), showWarnings = FALSE)

    ## Load datasets
    cellPopSpec = read_delim("cellPopSpec.txt", col_names = FALSE) %>%
        rename(g = 1, c = 2, s = 3, n = 4, biomass = 5, M = 6, pp = 7)
    totalPopSpec = read_delim("totalPopSpec.txt", col_names = FALSE) %>%
        rename(g = 1, s = 2, n = 3)
    traits = read_delim("../traits.txt", col_names = FALSE) %>%
        rename(M = 1, pp = 2) %>%
        add_column(.before = "M", s = 1:nrow(.))
    parameters = read_delim("../Parameters.txt", col_names = FALSE) %>%
        rename(name = 1, var = 2, value = 3)

    totalPop = totalPopSpec %>%
        group_by(g) %>%
        summarise(n = sum(n))

    ## Calculate total richness
    ## Only include species with a richness greater than 4
    totalRich = totalPopSpec %>%
        filter(n > 4) %>%
        group_by(g) %>%
        summarise(n = n_distinct(s))

    tClimate = parameters %>% filter(var == "climateT") %>% pull(value) %>% as.numeric()
    tMax = max(totalPopSpec$g)

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
    ggplot() + 
        geom_rect(aes(xmin = 0, xmax = tClimate, ymin = 0, ymax = Inf), fill = "lightblue", alpha = 0.5) +
        geom_rect(aes(xmin = tClimate, xmax = tMax, ymin = 0, ymax = Inf), fill = "red", alpha = 0.25) +
        geom_vline(xintercept = tClimate, linetype = "dashed", linewidth = 2) +
        geom_point(data = totalPop, aes(g, n), size = 2) +
        theme_classic() +
        labs(x = "Time", y = "Total Population") +
        theme(text = element_text(size = 30)) +
        scale_y_continuous(expand = c(0, 0)) +
        scale_x_continuous(expand = c(0, 0))

    ggsave(paste0(gpath, "../../../../Paper/Figures/climateChange/Seed_", x, "/totalPop_", x, ".pdf"), width = 15, height = 10)

    ## Plot total richness across the landscape over time
    ggplot() + 
        geom_rect(aes(xmin = 0, xmax = tClimate, ymin = 0, ymax = Inf), fill = "lightblue", alpha = 0.5) +
        geom_rect(aes(xmin = tClimate, xmax = tMax, ymin = 0, ymax = Inf), fill = "red", alpha = 0.25) +
        geom_vline(xintercept = tClimate, linetype = "dashed", linewidth = 2) +
        geom_point(data = totalRich, aes(g, n), size = 2) +
        theme_classic() +
        labs(x = "Time", y = "Species Richness") +
        theme(text = element_text(size = 30)) +
        scale_y_continuous(expand = c(0, 0)) +
        scale_x_continuous(expand = c(0, 0))

    ggsave(paste0(gpath, "../../../../Paper/Figures/climateChange/Seed_", x, "/totalRich_", x, ".pdf"), width = 15, height = 10)

    ## Plot SAD but with mass
    ggplot(filter(totalPopSpec, g == max(totalPopSpec$g)), aes(log10(M), log10(n), fill = as.factor(pp))) +
        geom_col(width = 0.01) +
        theme_classic() +
        theme(text = element_text(size = 30)) +
        labs(x = "Log10(Mass)", y = "Log10(Abundance", fill = NULL) +
        scale_y_continuous(expand = c(0, 0)) +
        scale_fill_viridis_d(labels = c("Consumer", "Primary Producer"))

    ggsave(paste0(gpath, "../../../../Paper/Figures/climateChange/Seed_", x, "/SAD_", x, ".pdf"), width = 15, height = 10)

    ## Test Damuth’s law (when logged the coefficient is the exponent)
    ggplot(totalPopSpec %>% filter (g == max(totalPopSpec$g) & n > 4), aes(log10(M), log10(n), col = as.factor(pp))) +
        geom_point(size = 7.5, alpha = 0.5) +
        geom_smooth(method = "lm", linewidth = 2) +
        theme_classic() +
        labs(x = "Log10(Body Mass)", y = "Log10(Abundance in MTaNa)", 
        col = NULL) +
        theme(text = element_text(size = 30)) +
        stat_poly_eq(use_label("eq"), size = 10, label.x = 0.9) +
        scale_colour_viridis_d(end = 0.7, labels = c("Consumer", "Primary Producer"))

    ggsave(paste0(gpath, "../../../../Paper/Figures/climateChange/Seed_", x, "/damuthLinear_", x, ".pdf"), width = 15, height = 10)


    ## Test damuth law over time
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
        group_split(g) %>%
        lapply(., checkDamuth) %>% 
        bind_rows()

    ## Plot damuth law over time
    ggplot() + 
        geom_rect(aes(xmin = 0, xmax = tClimate, ymin = -Inf, ymax = Inf), fill = "lightblue", alpha = 0.5) +
        geom_rect(aes(xmin = tClimate, xmax = tMax, ymin = -Inf, ymax = Inf), fill = "red", alpha = 0.25) +
        geom_vline(xintercept = tClimate, linetype = "dashed", linewidth = 2) +
        geom_line(data = damuth, aes(g, dam), linewidth = 1) +
        geom_line(linewidth = 1) +
        theme_classic() +
        geom_hline(yintercept = -0.75, linetype = "dashed", linewidth = 1) +
        theme(text = element_text(size = 30)) +
        labs(x = "Time", y = "Damuth's Exponent") +
        scale_y_continuous(expand = c(0, 0)) +
        scale_x_continuous(expand = c(0, 0))

    ggsave(paste0(gpath, "../../../../Paper/Figures/climateChange/Seed_", x, "/damuthTime_", x, ".pdf"), width = 15, height = 10)

    ## Plot the percentage of reproduction events that were within 1% of the maximum reproduction probability
    pMax = read_delim("pMax.txt", col_names = FALSE) %>%
        rename(g = 1, pMax = 2) %>%
        mutate(pProp = pMax/(max(g)/(max(g)/min(g))))

    ggplot() +
        geom_rect(aes(xmin = 0, xmax = tClimate, ymin = -Inf, ymax = Inf), fill = "lightblue", alpha = 0.5) +
        geom_rect(aes(xmin = tClimate, xmax = tMax, ymin = -Inf, ymax = Inf), fill = "red", alpha = 0.25) +
        geom_vline(xintercept = tClimate, linetype = "dashed", linewidth = 2) +
        geom_line(data = pMax, aes(g, pProp), linewidth = 1) +
        geom_line(linewidth = 1) +
        theme_classic() +
        theme(text = element_text(size = 30)) +
        labs(x = "Time", y = "Damuth's Exponent") +
        scale_y_continuous(expand = c(0, 0)) +
        scale_x_continuous(expand = c(0, 0))

}
