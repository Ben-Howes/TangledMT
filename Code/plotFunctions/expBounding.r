library(tidyverse)

vals = data.frame(H = seq(-10, 10, by = 0.1))
vals = vals %>% mutate(Poff = (exp(H))/(1 + exp(H)))

ggplot(vals, aes(H, Poff)) + 
    geom_line(size = 3) +
    theme_classic() +
    labs(x = expression("H/M (kg s"^"-1"~"kg"^"-1)"), y = "Poff") +
    theme(text = element_text(size = 30))

ggsave(filename = "expFigureHM.png", width = 15, height = 8)
