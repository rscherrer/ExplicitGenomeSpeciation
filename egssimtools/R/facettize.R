#' Flexible plot facetting
#'
#' Customizable way to split a ggplot into multiple facets depending on facetting variables.
#'
#' @param p A ggplot
#' @param data Data frame used to produce the plot, here used to read the labels
#' @param facet_rows Optional parameter name(s) to facet the plot by rows
#' @param facet_cols Optional parameter name(s) to facet the plot by columns
#' @param facet_wrapped Whether to automatically fill the facets by row with one-dimensional array of parameter combinations. If TRUE, combines the parameters in facet_rows and facet_cols into a single array. Ignored if none of those is specified.
#' @param label_facets Whether to add custom labels to the facets
#' @param facet_prefixes Optional prefixes to add to the facet labels of each parameter. Must have one value per facetting parameter. Ignored if label_facets is FALSE. If not specified and label_facets is TRUE, the names of the parameters are used as prefixes.
#' @param sep Optional separator to use when adding prefixes. Defaults to " = ".
#'
#' @return A facetted ggplot
#'
#' @export

facettize <- function(
  p,
  data,
  facet_rows = NULL,
  facet_cols = NULL,
  facet_wrapped = NULL,
  label_facets = FALSE,
  facet_prefixes = NULL,
  sep = " = "
) {

  library(tidyverse)

  facets <- c(facet_rows, facet_cols)
  if (is.null(facets)) return (p)

  # Setup rows and columns
  if (!is.null(facet_rows)) lhs <- paste(facet_rows, collapse = " + ") else lhs <- "."
  if (!is.null(facet_cols)) rhs <- paste(facet_cols, collapse = " + ") else rhs <- "."

  # Choose what type of facetting to use
  if (facet_wrapped) thisfacet <- facet_wrap else thisfacet <- facet_grid

  # Add custom labels if needed
  if (is.null(facet_prefixes)) facet_prefixes <- facets
  labels <- make_facet_labels(data, facets, facet_prefixes, sep = sep, no_use = !label_facets, simplify = FALSE)

  # Split the plot into facets
  p + thisfacet(formula(paste(lhs, "~", rhs)), labeller = do.call("labeller", labels))

}
