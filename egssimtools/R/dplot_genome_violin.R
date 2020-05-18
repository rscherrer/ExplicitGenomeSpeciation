#' Diagnostic violin plot
#'
#' Quick-plots a variable across the genome, possibly split between categories
#'
#' @param root Path to the simulation folder
#' @param y Name of the variable to plot
#' @param x Optional grouping variable
#' @param t Time point at which to show the scan. If unspecified, last time point.
#'
#' @return A ggplot
#'
#' @export

dplot_genome_violin <- function(root, y, x = NULL, t = NULL) {

  library(ggplot2)

  time <- read_data(root, "time")
  X <- read_data(root, y)
  nloci <- nrow(X) / nrow(time)

  data <- read_data(root, c("time", y), dupl = c(nloci, 1), architecture = TRUE)

  if (is.null(t)) t <- last(data$time)

  rmaxis <- FALSE
  if (is.null(x)) {
    x.. <- "1"
    x <- "x.."
    rmaxis <- TRUE
  }

  data <- data %>% filter(time == t)
  p <- ggplot(data, aes(x = factor(get(x)), y = get(y))) +
    geom_violin() +
    theme_bw() +
    labs(x = x, y = y)

  if (rmaxis) p <- p + theme(
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_blank()
  )

  return(p)

}
