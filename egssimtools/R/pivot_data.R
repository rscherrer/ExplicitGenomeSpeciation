#' Pivot data to long format
#'
#' Can be useful when plotting with facets.
#'
#' @param data A data frame
#' @param variables The columns to pass as `cols` in `tidyr::pivot_longer`
#' @param newnames Optional character vector of new names for `variables`
#'
#' @return A tibble
#'
#' @examples
#'
#' \dontrun{
#'
#' # Location of the simulation data
#' root <- "data/example_1"
#'
#' # Read speciation metrics in the long format
#' data <- read_sim(root, c("EI", "RI", "SI"))
#' pivot_data(data, c("EI", "RI", "SI"))
#'
#' # Read individual trait values in the long format
#' data <- read_pop(root, "individual_trait", by = 3)
#' cols <- paste0("individual_trait", 1:3)
#' pivot_data(data, cols, newnames = 0:2)
#'
#' }
#'
#' @export

pivot_data <- function(data, variables, newnames = NULL) {

  data <- tidyr::pivot_longer(
    data,
    cols = tidyselect::all_of(variables),
    names_to = "variable",
    values_to = "value"
  )

  if (!is.null(newnames)) {

    names(variables) <- newnames
    data <- data %>%
      dplyr::mutate(variable = factor(variable)) %>%
      dplyr::mutate(variable = forcats::fct_recode(variable, !!!variables))

  }

  return(data)

}
