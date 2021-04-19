#' Plot genomic heatmap through time
#'
#' Heatmap of a variable through time across the genome.
#'
#' @param data A locus-wise data frame with columns "time", "locus" and `z` at least
#' @param z Name of the variable to plot
#'
#' @return A ggplot
#'
#' @export

plot_genome_heatmap <- function(

  data,
  z

) {

  # Set up
  p <- ggplot2::ggplot(data, ggplot2::aes(time, y = locus, fill = get(z))) +
    geom_tile() +
    labs(fill = z)

  return(p)

}
