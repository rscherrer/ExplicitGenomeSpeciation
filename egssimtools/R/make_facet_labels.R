#' Make facet labels
#'
#' Make facet labels of the form "variable = value".
#'
#' @param data The dataset
#' @param columns The variables to make facet labels for
#' @param prefixes Prefixe to add to each facet variable. Can be a named vector where the names are the corresponding facetting variables. Must be the same length as columns.
#' @param no_use Whether to return labels without prefix. Useful only in specific cases e.g. use inside other functions.
#' @param sep Separator between prefix and label. Defaults to " = ".
#' @param simplify If one facetting variable is supplied, whether to return a vector rather than a list of one vector
#'
#' @return A named list of named vectors of labels for each facetting variable. The output is ready for use as the "labeller" argument of the facet_wrap and facet_grid functions in ggplot by using labeller = do.call("labeller", labels), where labels is the output of this function. Returns a vector of named labels if only one facetting variable is requested and simplify is TRUE. In this case use labeller = labeller(ecosel = labels), assuming that labels is the function output and ecosel is the name of the facetting variable.
#'
#' @export

make_facet_labels <- function(data, columns, prefixes = NULL, no_use = FALSE, sep = " = ", simplify = TRUE) {

  library(assertthat)

  assert_that(all(columns %in% colnames(data)))

  # Read levels for each facetting column
  labels <- lapply(columns, function(column) {
    column <- unlist(data[, column])
    assert_that(is.factor(column))
    labels <- levels(column)
    names(labels) <- labels
    return (labels)
  })
  names(labels) <- columns

  # Return as is if needed (specific cases like use in other functions)
  if (no_use) if (simplify & length(columns) == 1) return (labels[[1]]) else return (labels)

  # Set prefixes if needed
  if (!is.null(prefixes)) {

    if (length(prefixes) != length(columns)) stop("Please provide a prefix for each facetting column")
    if (is.null(names(prefixes))) columns <- prefixes else columns <- prefixes[columns]

  }

  # Add prefixes to the labels
  labels <- mapply(function(labels, column) gsub("^", paste0(column, sep), labels), labels, columns, SIMPLIFY = FALSE)

  if (simplify & length(columns) == 1) labels[[1]] else labels

}
