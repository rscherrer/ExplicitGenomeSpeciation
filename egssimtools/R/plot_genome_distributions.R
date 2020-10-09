#' Plot distributions across loci per trait
#'
#' @param data A locus-wise dataset.
#' @param y Name of the variable to plot.
#' @param t Optional time point to filter. Defaults to the last one if `timed` is TRUE.
#' @param timed Whether there are multiple time points in the dataset (the column must be named "time").
#' @param jitter Logical. Whether to show jittered points.
#' @param width Width of the jitter.
#'
#' @return A ggplot
#'
#' @export

plot_genome_distributions <- function(

  data,
  y,
  t = NULL,
  timed = TRUE,
  jitter = TRUE,
  width = 0.2

) {

  if (timed) {

    if  (is.null(t)) t <- dplyr::last(unique(data$time))
    data <- data %>% dplyr::filter(time == t)

  }

  p <- ggplot2::ggplot(data, ggplot2::aes(x = trait, y = get(y), color = trait)) +
    ggplot2::geom_violin(draw_quantiles = 0.5)

  if (jitter) p <- p + ggplot2::geom_jitter(width = 0.2)

  return(p)

}
