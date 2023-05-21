#################################################
## Analyse TaNa output where it only has
## trophic interactions
#################################################

library(tidyverse)

gpath = "/home/ben/Documents/TangledMT/Results/TNM_Output/Seed_2/Results/"
setwd(gpath)

totalPop = read_delim("totalPop.txt", col_names = FALSE) %>%
    rename(g = 1, n = 2)


ggplot(totalPop, aes(g, n)) + geom_line() +
    theme_classic()

totalPopSpec = read_delim("totalPopSpec.txt", col_names = FALSE) %>%
    rename(g = 1, s = 2, n = 3)

ggplot(filter(totalPopSpec, n > 100), aes(g, n, col = as.factor(s))) + geom_line() +
    theme_classic()

ggplot(filter(totalPopSpec, g == max(totalPopSpec$g)), aes(fct_rev(fct_reorder(as.factor(s), n)), n)) +
    geom_col() +
    theme_classic()
