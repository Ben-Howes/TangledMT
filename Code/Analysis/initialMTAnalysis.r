#################################################
## Analyse TaNa output where it only has
## trophic interactions
#################################################

library(tidyverse)

gpath = "/home/ben/Documents/TangledMT/Results/TNM_Output/Seed_5/Results/"
setwd(gpath)

## Load datasets

totalPop = read_delim("totalPop.txt", col_names = FALSE) %>%
    rename(g = 1, n = 2)
cellPopSpec = read_delim("cellPopSpec.txt", col_names = FALSE) %>%
    rename(g = 1, c = 2, s = 3, n = 4, H = 5)
totalPopSpec = read_delim("totalPopSpec.txt", col_names = FALSE) %>%
    rename(g = 1, s = 2, n = 3)
traits = read_delim("../Traits.txt", col_names = FALSE) %>%
    rename(M = 1, pp = 2) %>% dplyr::select(M) %>%
    add_column(.before = "M", s = 1:nrow(.))

## Join trait data with totalpopspec
totalPopSpec = totalPopSpec %>% left_join(traits)

ggplot(totalPop, aes(g, n)) + 
    geom_line() +
    theme_classic() +
    labs(x = "Generation", y = "Total Population") +
    theme(text = element_text(size = 30))

ggplot(filter(totalPopSpec, g == max(totalPopSpec$g)), aes(fct_rev(fct_reorder(as.factor(M), n)), n)) +
    geom_col() +
    theme_classic() +
    theme(text = element_text(size = 30),
    axis.text.x=element_blank()) +
    labs(x = "Body Mass", y = "Abundance") +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 2000)) +
    geom_text(aes(label = round(M/1000)), position = position_dodge(width=0.9), vjust=-0.25, size = 5)

### Change in mean mass over time, calculated as (mass*abundance)/abundance
avgMass = totalPopSpec %>% 
    mutate(totalM = n*M) %>%
    group_by(g) %>%
    summarise(avgM = sum(totalM)/sum(n))

ggplot(avgMass, aes(g, avgM/1000)) + 
    geom_line() +
    theme_classic() +
    labs(x = "Time", y = "Average Mass (kg)") +
    theme(text = element_text(size = 30))


### Explore cellPopSpec
## Join cellPopSpec with traits

cellPopSpec = cellPopSpec %>% left_join(traits)
test = cellPopSpec %>% filter(g == max(cellPopSpec$g) & c == 314)
test
