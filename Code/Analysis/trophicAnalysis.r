#################################################
## Analyse TaNa output where it only has
## trophic interactions
#################################################

library(tidyverse)

gpath = "/home/ben/Documents/TangledMT/Results/TNM_Output/Seed_2/Results/"
setwd(gpath)

totalPop = read_delim("totalPop.txt", col_names = FALSE) %>%
    rename(g = 1, n = 2)


ggplot(totalPop, aes(g, n)) + 
    geom_line() +
    theme_classic() +
    labs(x = "Generation", y = "Total Population") +
    theme(text = element_text(size = 30))

ggsave(filename = "../../../../Paper/Figures/abundance.png", width = 10, height = 6)

totalPopSpec = read_delim("totalPopSpec.txt", col_names = FALSE) %>%
    rename(g = 1, s = 2, n = 3)

ggplot(filter(totalPopSpec, g == max(totalPopSpec$g)), aes(fct_rev(fct_reorder(as.factor(s), n)), n)) +
    geom_col() +
    theme_classic() +
    theme(text = element_text(size = 30),
    axis.text.x = element_blank()) +
    labs(x = "Species", y = "Abundance") +
    scale_y_continuous(expand = c(0, 0))

ggsave(filename = "../../../../Paper/Figures/SADs.png", width = 10, height = 6)
