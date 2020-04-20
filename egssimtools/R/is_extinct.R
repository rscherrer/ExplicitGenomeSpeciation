#' Did a simulation go extinct?
#'
#' @param folder Path to the folder
#'
#' @export

is_extinct <- function(folder) {

  # Look for the word "extinct" in the SLURM output file

  i <- grep("slurm.*out", list.files(folder))
  if (length(i) == 0) return (FALSE)
  if (length(i) > 1) warning(paste0("More than one SLURM output file found in ", folder, ", taking the last one."))
  i <- i[length(i)]
  filename <- list.files(folder, full.names = TRUE)[i]
  lines <- read.delim(filename, header = FALSE)[, 1]
  length(grep("extinct", lines)) > 0

}
