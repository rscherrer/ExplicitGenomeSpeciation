# Use this script to set-up parameter files for many combinations
# and move the files into their respective folders

rm(list = ls())

library(egssimtools)
library(tidyr)

# Setup parameter combinations
hsymmetry <- c(0, 1)
scaleI <- c('0 0 0', '1 1 1')
seed <- c(24, 42, 53, 12, 66)

cmb <- expand_grid(hsymmetry, scaleI, seed)

# Add extra parameters
cmb <- cmb %>%
  mutate(
    ecosel = ifelse(hsymmetry == 0, 1, 1.6),
    scaleA = ifelse(scaleI == "0 0 0", "1 1 1", "0 0 0"),
    tend = 20000,
    tsave = 100
  )

# Make parameter files and save them into respective folders
cmb <- cmb %>% split(f = factor(seq(nrow(cmb))))
folders <- paste0("sim", 1:20)
folders %>% map(dir.create)
filename <- "%s/parameters.txt"
template <- "parameters.txt"

map2(cmb, folders, ~ set_param_file(.x, template, saveto = sprintf(filename, .y)))

# Write a whattosave into each simulation folder
whattosave <- c(
  "time",
  "EI",
  "RI",
  "SI",
  "population_size",
  "individual_trait",
  "Fst",
  "Gst",
  "Qst",
  "Cst",
  "varG",
  "varA",
  "varN",
  "varT",
  "genome_Fst",
  "genome_Gst",
  "genome_Qst",
  "genome_Cst",
  "genome_varG",
  "genome_varA",
  "genome_varN",
  "genome_alpha",
  "genome_meang",
  "genome_freq",
  "network_corgen",
  "network_corbreed",
  "network_corfreq",
  "network_avgi",
  "network_avgj"
)

map(
  sprintf("%s/whattosave.txt", folders),
  function(fname) {
    f <- file(fname)
    writeLines(whattosave, f)
    close(f)
  }
)
