#' Genome scan plot
#'
#' Flexible function to plot a variable across loci along the genome.
#' The data can be plotted for a given time point in a simulation,
#' or for no specific time point if the data are not timed (e.g. summary data).
#' The plot can be animated through time points. Edges of the gene network
#' can be displayed onto the layer of points representing the loci.
#' Supports additional aesthetics mapping and customization.
#'
#' @param data A data frame containing information on a per locus (and optionally per time point) basis
#' @param mode Either of "timepoint" to filter a specific time point, "summary" if there are no multiple time points, and "animation" if an animation through time points should be made
#' @param t Integer. The time point to filter
#' @param x String. Either of "location" to show loci locations along the genome as the horizontal axis, or "locus" to show loci indices
#' @param y String. Name of the variable on the vertical axis
#' @param show_generation Logical. Whether to show the time point as a title
#' @param style String. The `geom`s to use. Either of "points", "bars" or "lollipop".
#' @param color Name of the variables to map to those aesthetics.
#' @param chrom_pars List of graphical parameters for the chromosomes, to be passed to `add_chromosomes()`. Includes `y`, `offset`, `size` and `color`.
#'
#' @return A ggplot, or an animation from the `gganimate` package
#'
#' @note This is a simplified version of `plot_genome_scatter()`
#'
#' @export

plot_genome_scan <- function(

  data,
  mode = "timepoint",
  t = NULL,
  x = "location",
  y = "genome_Fst",
  show_generation = FALSE,
  style = "lollipop",
  color = NULL,
  show_chroms = FALSE,
  chrom_pars = list(y = -0.5, offset = -0.01, size = 2, color = "black")

) {

  if (mode == "summary") show_generation <- FALSE

  # Filter a single time point if needed
  if (mode == "timepoint") {

    if (is.null(t)) t <- dplyr::last(unique(data$time))
    data <- dplyr::filter(data, time == t)

  }

  # Set up
  p <- ggplot2::ggplot(data)

  # Prepare aesthetic mapping
  aes_string <- "x = get(x), y = get(y)"
  if (!is.null(color)) aes_string <- paste0(aes_string, ", color = get(color)")

  # Add geoms (bars, points or both = lollipops)
  if (style %in% c("bars", "lollipop")) {

    aes_string_segments <- paste0(aes_string, ", xend = get(x)")
    aes_string_segments <- paste0("ggplot2::aes(", aes_string_segments, ")")

    p <- p +
      ggplot2::geom_segment(
        mapping = eval(rlang::parse_expr(aes_string_segments)),
        yend = 0
      )
  }

  if (style %in% c("points", "lollipop")) {

    aes_string <- paste0("ggplot2::aes(", aes_string, ")")
    p <- p + ggplot2::geom_point(mapping = eval(rlang::parse_expr(aes_string)))

  }

  # Custom axis labels
  p <- p + xlab(x) + ylab(y)

  # Show the generation as title
  if (show_generation) {
    title <- paste(
      "Generation",
      ifelse(mode == "timepoint", t, "{round(frame_time, -2)}")
    )
    p <- p + ggplot2::labs(title = title)
  }

  # Show chromosomes as horizontal bars
  if (show_chroms) {
    p <- add_chromosomes(
      p, y = chrom_pars$y, offset = chrom_pars$offset, size = chrom_pars$size,
      color = chrom_pars$color, unit = x
    )
  }

  # Animate if needed
  if (mode == "animation") p <- p + gganimate::transition_time(time)

  return(p)

}
