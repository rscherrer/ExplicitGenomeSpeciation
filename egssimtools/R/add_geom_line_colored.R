#' Customize line colors
#'
#' Easy way to add lines to a ggplot and color them with custom colors and variables
#'
#' @param p A ggplot
#' @param color_by The variable to use for coloration
#' @param color_by_numeric Whether to convert the color variable into a numeric
#' @param colors Optional colors of the lines. If color_by is NULL, the color that all the lines must have. If color_by is defined and color_by_numeric is TRUE, the lower and the upper end of the color gradient. If color_by_numeric is FALSE, a set of colors for each level of the color factor. If not defined, default ggplot color(s) will be used.
#'
#' @param return A ggplot with colored lines added
#'
#' @export

add_geom_line_colored <- function(p, color_by = NULL, color_by_numeric = FALSE, colors = NULL) {

  library(tidyverse)

  if (is.null(color_by)) return (p + geom_line(color = ifelse(is.null(colors), "black", colors[1])))

  # Color according to a parameter
  if (color_by_numeric) thisconvert <- function(x) as.numeric(as.character(x)) else thisconvert <- factor
  p <- p + geom_line(aes(color = thisconvert(get(color_by)))) + labs(color = color_by)

  if (!is.null(colors)) {

    # Add a custom color scale if provided
    if (color_by_numeric) {

      # Either a color gradient
      if (length(colors) != 2) stop("Please provide a lower and upper end for color gradient")
      p <- p + scale_color_gradient(low = colors[1], high = colors[2])

    } else {

      # Or a custom set of colors
      p <- p + scale_color_manual(values = colors)

    }
  }

  return (p)

}
