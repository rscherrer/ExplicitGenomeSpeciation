rm(list = ls())

library(egssimtools)
library(tidyverse)

root <- "/media/raphael/bigass/simulations/EGS/genomes6/batch4/sim_ecosel_1.4_hsymmetry_1_scaleA_0_0_0_scaleI_1_1_1_archfile_arch1.txt_r2/"

data <- read_indiv_loci(root, t = 19900, extra = "individual_ecotype")

data <- data %>% filter(locus == 1)
ggplot(data, aes(x = factor(allcount), y = value, color = factor(ecotype))) +
  geom_violin()
