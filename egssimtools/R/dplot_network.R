#' Plot gene network
#'
#' @param root Path to the simulation folder
#' @param y Name of the variable to map as node size aesthetics
#' @param t What time point? Defaults to last time point
#'
#' @return A ggplot
#'
#' @export

dplot_network <- function(root, y, t = NULL) {

  library(tidyverse)
  library(ggraph)
  library(tidygraph)

  time <- read_data(root, "time")
  X <- read_data(root, y)
  nloci <- nrow(X) / nrow(time)

  data <- read_data(root, c("time", y), dupl = c(nloci, 1), architecture = TRUE)
  if (is.null(t)) t <- last(data$time)
  data <- data %>% filter(time == t)

  arch <- read_network_architecture(root)
  arch <- arch %>% activate(nodes) %>% right_join(data)
  arch <- arch %>% activate(nodes) %>% mutate(trait = factor(trait))
  arch <- arch %>% activate(edges) %>% mutate(trait = factor(trait))

  colors <- c("forestgreen", "goldenrod", "grey")
  p <- ggraph(arch, layout = "graphopt", charge = 0.1, mass = 30, niter = 100000) +
    geom_edge_link(aes(color = factor(trait)), width = 2, alpha = 0.6) +
    geom_node_point(aes_string(fill = "trait", size = y), shape = 21) +
    scale_edge_color_manual(values = colors) +
    scale_fill_manual(values = colors) +
    scale_size_continuous(range = c(2, 7)) +
    scale_alpha(range = c(0.6, 1)) +
    labs(size = y, fill = "Trait") +
    theme_void() +
    guides(edge_color = FALSE)

  p

}
