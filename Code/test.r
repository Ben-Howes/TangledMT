library(tidyverse)

path = "/home/ben/Documents/TangledMT/Results/TNM_Output/Seed_2/Results/"
cellPopSpec = read_delim(paste0(path, "cellPopSpec.txt"), col_names = F) %>%
rename(g = 1, c = 2, s = 3, n = 4)
totalPop = read_delim(paste0(path, "totalPop.txt"), col_names = F) %>%
rename(g = 1, n = 2)
cellPop = read_delim(paste0(path, "cellPop.txt"), col_names = F) %>%
rename(g = 1, c = 2, n = 3)
totalPopSpec = read_delim(paste0(path, "totalPopSpec.txt"), col_names = F) %>%
rename(g = 1, s = 2, n = 3)
totalRich = read_delim(paste0(path, "totalRich.txt"), col_names = F) %>%
rename(g = 1, n = 2)

cellPopSpec %>% filter(g == 1000) %>% summarise(N = sum(n))
totalPop %>% filter(g == 1000)
cellPopSpec %>% filter(g == 1000) %>% group_by(c) %>% summarise(N = sum(n))
cellPop %>% filter(g == 1000) %>% arrange(c)
cellPopSpec %>% filter(g == 1000) %>% group_by(s) %>% summarise(N = sum(n))
totalPopSpec %>% filter(g == 1000) %>% arrange(s)
cellPopSpec %>% filter(g == 1000) %>% distinct(s) %>% nrow()
totalRich %>% filter(g == 1000)
