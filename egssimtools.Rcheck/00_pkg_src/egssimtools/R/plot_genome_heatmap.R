#' Genome heatmap
#'
#' Plots a heatmap of a variable across the genome, through time
#'
#' @param folder The simulation folder
#' @param variable What variable to plot
#'
#' @return A ggplot
#'
#' @export

plot_genome_heatmap <- function(folder, variable) {

  library(tidyverse)

  data <- read_loci(folder, variable)
  data <- do.call("cbind", data) %>%
    as.data.frame() %>%
    gather("time", "variable") %>%
    mutate(time = as.numeric(time)) %>%
    group_by(time) %>%
    mutate(locus = seq(n()))

  ggplot(data, aes(x = time, y = locus, fill = variable)) +
    geom_tile() +
    labs(fill = variable) +
    scale_fill_continuous(type = "viridis") +
    theme_bw()

}
