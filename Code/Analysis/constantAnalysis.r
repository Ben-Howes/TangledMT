############################################
## Script to analyse outputs of MTaNa
## using different constants
############################################

library(tidyverse)
library(ggpmisc) ## stat_poly
library(lemon) ##facet_rep_wrap

gpath = "/home/ben/Documents/TangledMT/Results/TNM_Output/"
setwd(gpath)

## We want to analyse simple metrics between models run with different constants
## I plan to test differences in:
## 1) Average total abundance, 2) Variation in total abundance
## 3) Average body mass, 4) Variation in body mass
## 5) Average dalmuth rule constant (and SE)

## Get paths to files from each run
paths = list.files(gpath)

## Read in parameters, seed, and cellPopInd for each run - join cellPopInd with parameters

getCellPopInd = function(path) {
    params = read_table(paste0(gpath, path, "/Parameters.txt"), col_names = F) %>%
        dplyr::select(X2, X3) %>% rename("constant" = 1, "value" = 2) %>%
        pivot_wider(names_from = "constant", values_from = "value")
    cellPopInd = read_table(paste0(gpath, path, "/Results/cellPopInd.txt"), col_names = F) %>%
        rename("g" = 1, "c" = 2, "s" = 3, "m" = 4, "pp" = 5, "disp" = 6) %>%
        mutate(.before = "g", seed = path)
    cellPopInd = cbind(cellPopInd, params)

    return(cellPopInd)
}

cellPopInd = lapply(paths[1:3], getCellPopInd) %>% bind_rows()

totalPop = cellPopInd %>% group_by(seed, g, I0) %>% 
    summarise(genPop = n()) %>% 
    group_by(seed, I0) %>%
    summarise(avgTotalPop = mean(genPop), sdTotalPop = sd(genPop))
