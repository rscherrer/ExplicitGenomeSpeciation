#' Smoothen lines using LOESS
#'
#' @param data A data frame
#' @param x Name of the horizontal axis variable
#' @param y Name of the vertical axis variable
#' @param span Span parameter for the smoothing
#' @param line Name of the variable identifying the different lines. Leave
#' unspecified if there is only one line to smoothen.
#'
#' @return A data frame
#'
#' @examples
#'
#' root <- system.file("extdata", "example_1", package = "egssimtools")
#' data <- read_loci(root, "genome_Fst", architecture = TRUE)
#' smoothen_data(data, "time", "genome_Fst", line = "locus")
#'
#' @export

smoothen_data <- function(data, x, y, span = 0.2, line = NULL) {

  smoothen <- function(data) {
    loess(
      as.formula(paste(y, "~", x)), degree = 1, span = span, data = data
    )$fitted
  }

  if (is.null(line)) {
    data <- data %>% dplyr::mutate(linecol = 1)
    line <- "linecol"
  }

  data <- data %>%
    dplyr::group_by_at(line) %>%
    tidyr::nest() %>%
    dplyr::mutate(smooth = purrr::map(data, smoothen)) %>%
    tidyr::unnest(cols = c(data, smooth)) %>%
    dplyr::ungroup()

  if ("linecol" %in% colnames(data)) data <- data %>% dplyr::select(-linecol)

  return(data)

}
