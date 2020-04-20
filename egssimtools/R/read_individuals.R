#' Read individual data
#'
#' @param filename Path to the folder
#' @param variable What to read
#' @param nvalues How many values per individual? 3 for traits, 1 for habitat and ecotype
#'
#' @export

read_individuals <- function(folder, variable, nvalues = 3) {

  data <- read_data(folder, variable)

  t <- read_time(folder)
  n <- read_population_size(folder)

  # Assign each data value to a time point
  t <- rep(mrep(t, n), each = nvalues)

  # And to an individual within time points
  n <- lapply(n, function(n) rep(seq_len(n), each = nvalues))

  # Split the data by individual and by time point
  mapply(function(data, n) split(data, n), split(data, t), n, SIMPLIFY = FALSE)

}
