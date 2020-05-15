#' Read time only if available
#'
#' Read the time.dat file only if it is present. If not, return a sequence of time point indices.
#'
#' @param filename Path to the folder
#' @param ntimes Number of time points if the time file is not present
#'
#' @export

read_time_if <- function(folder, ntimes = NULL) {

  if ("time.dat" %in% list.files(folder)) return (read_time(folder))
  if (is.null(ntimes)) stop("Please provide a number of time points")
  return (seq_len(ntimes))

}
