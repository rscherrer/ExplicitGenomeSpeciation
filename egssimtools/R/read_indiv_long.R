#' Read individual data through time (long format)
#'
#' Wrapper around `read_data` to read individual-wise data in the long format
#' (useful for facetting plots).
#'
#' @param folder Path to the folder
#' @param variables What variables to read (`time` is included by default)
#' @param by See `?read_indiv`.
#' @param to_long What variable (within `variables`) to use to pivot the data
#' in the long format? Defaults to the first element in `variables`.
#' See `tidyr::pivot_longer` for details.
#' @param parnames,combine,as_numeric,parfile Parameters for `read_data`
#' @param names_to Optional name for the new column containing variable names.
#' See `?tidyr::pivot_longer` for details.
#' @param replace_with Optional string to replace `to_long` in the column that
#' will contain the variable names (e.g. change from "individual_trait1" to
#' "trait 1"). No replacement if NULL.
#' @param cpp_numbering Whether to keep the C++ numbering style in variable
#' names (starting from 0)
#'
#' @details The files `time.dat` and `population_size.dat` must be present.
#' Multiple variables can be read, but one has to be chosen for pivoting to the
#' long format.
#'
#' @return A data frame
#'
#' @examples
#'
#' root <- system.file("extdata", "example_1", package = "egssimtools")
#' read_indiv_long(root, "individual_trait", by = 3, replace_with = "trait ")
#'
#' @export

read_indiv_long <- function(
  folder,
  variables,
  by = 3,
  to_long = variables[1],
  parnames = NULL,
  combine = FALSE,
  as_numeric = NULL,
  parfile = "paramlog.txt",
  names_to = "variable",
  replace_with = "",
  cpp_numbering = TRUE
) {

  if (!to_long %in% variables) stop("to_long must be one of variables")

  data <- read_data(
    folder,
    c("time", variables),
    by = c(1, by),
    dupl = list("population_size", 1),
    parnames = parnames,
    combine = combine,
    as_numeric = as_numeric,
    parfile = parfile
  )

  ycols <- colnames(data)[grep(to_long, colnames(data))]

  data <- data %>%
    tidyr::pivot_longer(
      cols = all_of(ycols),
      names_to = names_to,
      values_to = "value"
    )

  if (cpp_numbering) {
    data <- data %>%
      dplyr::mutate(
        extra = as.character(
          as.numeric(stringr::str_replace(get(names_to), to_long, "")) - 1
        )
      ) %>%
      dplyr::mutate_at(
        names_to,
        ~ stringr::str_replace(.x, "[0-9]", extra)
      ) %>%
      dplyr::select(-extra)
  }

  if (!is.null(replace_with)) {
    data <- data %>%
      dplyr::mutate_at(
        names_to,
        ~ stringr::str_replace(.x, to_long, replace_with)
      )
  }

  return(data)

}
