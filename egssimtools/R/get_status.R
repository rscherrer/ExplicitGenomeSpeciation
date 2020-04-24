#' SLURM status
#'
#' Read the exit status of the simulation from the SLURM output file
#'
#' @param folder Path to the folder
#'
#' @export

get_status <- function(folder) {

  i <- grep("slurm.*out", list.files(folder))
  if (length(i) == 0) return (NULL)
  if (length(i) > 1) warning(paste0("More than one SLURM output file found in ", folder, ", taking the last one."))
  i <- i[length(i)]
  filename <- list.files(folder, full.names = TRUE)[i]
  lines <- read.delim(filename, header = FALSE)[, 1]
  gsub("^.*: ", "", as.character(lines[grep("^State.*:", lines)]))

}
