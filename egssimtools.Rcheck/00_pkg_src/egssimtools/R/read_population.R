#' Read population-wide data
#'
#' @param filename Path to the folder
#' @param variable What to read
#' @param nvalues Number of values per time point
#' @param add_time Whether to try to append a time vector to the result
#'
#' @export

read_population <- function(folder, variable, nvalues = 1, add_time = FALSE) {

  data <- read_data(folder, variable)

  if (nvalues == 1) return(data)

  ntimes <- length(data) / nvalues
  t <- rep(seq(ntimes), each = nvalues)
  data <- split(data, t) %>%
    do.call("rbind", .) %>%
    data.frame %>%
    rename_all(str_replace, "X", variable)

  if (add_time) data <- data %>% mutate(time = read_time_if(folder, ntimes))

  return (data)

}
