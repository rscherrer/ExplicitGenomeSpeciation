rm(list = ls())

library(egssimtools)
library(tidyverse)

data <- lapply(sprintf("/media/raphael/bigass/simulations/EGS/EGS_sim%s", seq(6)), function(root) {
  collect_simulations(
    root,
    c("EI", "SI", "RI"),
    parnames = c("ecosel", "hsymmetry", "dispersal", "mutation", "scaleA", "scaleI"),
    to_numeric = c("ecosel", "hsymmetry", "dispersal", "mutation"),
    as_address = TRUE
  )
})
data <- do.call("rbind", data)

saveRDS(data, "data/population_wide_data.rds")
