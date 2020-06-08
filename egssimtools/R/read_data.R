#' Read simulation data
#'
#' Function to read simulation data from binary files into a data frame. Can read multiple data files and bind them together.
#'
#' @param folder Path to the folder
#' @param variables What variables to read
#' @param by A list. For each variable, length of the chunks to split it by. Recycled if its length is one.
#' @param dupl A list. For each variable, how many times to duplicate each row. Provide an integer to duplicate every row the same number of times, a vector of integers to duplicate each row a specific number of times, or a character string to read a vector of number of times from a .dat file. This argument is recycled if its length is one.
#' @param parnames Character vector of parameter names to append to the simulation data
#' @param combine,as_numeric Parameters for `read_parameters`
#' @param architecture Whether to read and append locus-wise genetic architecture. Make sure then that you are reading `variables` that are locus-specific variables in long-format (one column per variable).
#' @param archfile Name of the architecture file, if needed
#' @param parfile Name of the parameter file, if needed
#'
#' @return A data frame
#'
#' @note Do not provide the extension of the files. It is assumed to be `.dat`.
#'
#' @examples
#'
#' root <- system.file("extdata", "example_1", package = "egssimtools")
#' read_data(root, c("time", "genome_Fst"), by = c(1, 300))
#'
#' @export

read_data <- function(
  folder,
  variables,
  by = 1,
  dupl = 1,
  parnames = NULL,
  combine = FALSE,
  as_numeric = NULL,
  architecture = FALSE,
  archfile = "architecture.txt",
  parfile = "paramlog.txt"
) {

  data <- list(variables, by, dupl) %>%
    purrr::pmap_dfc(function(variable, by, dupl) {

      data <- read_binary(paste0(folder, "/", variable, ".dat")) %>%
        split(rep(seq(length(.) / by), each = by)) %>%
        do.call("rbind", .) %>%
        data.frame
      if (ncol(data) > 1) {
        colnames <- paste0(variable, seq(ncol(data)))
      } else {
        colnames <- variable
      }
      if (is.character(dupl)) {
        dupl <- read_binary(paste0(folder, "/", dupl, ".dat"))
      }
      if (length(dupl) == 1) dupl <- rep(dupl, nrow(data))
      data <- data[mrep(seq(nrow(data)), n = dupl), ] %>% data.frame
      data <- data %>% rename_str(colnames)

    })

  if (!is.null(parnames)) {

    pars <- read_parameters(
      folder, parnames, combine = combine, flatten = TRUE,
      as_numeric = as_numeric, filename = parfile
    )
    pars <- pars %>% purrr::map_dfc(rep, nrow(data))
    data <- list(data, pars) %>% dplyr::bind_cols

  }

  if (architecture) {

    arch <- read_genome_architecture(folder, filename = archfile)
    ntimes <- nrow(data) / nrow(arch)
    arch <- purrr::map_dfr(seq(ntimes), ~ arch)
    data <- cbind(data, arch)

  }

  return (data)

}
