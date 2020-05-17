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
#'
#' @return A data frame
#'
#' @note Do not provide the extension of the files. It is assumed to be .dat.
#'
#' @export

read_data <- function(
  folder, variables, by = 1, dupl = 1, parnames = NULL, combine = FALSE,
  as_numeric = NULL
) {

  library(tidyverse)

  data <- list(variables, by, dupl) %>%
    pmap_dfc(function(variable, by, dupl) {

      data <- read_binary(paste0(folder, "/", variable, ".dat")) %>%
        split(rep(seq(length(.) / by), each = by)) %>%
        do.call("rbind", .) %>%
        data.frame
      if (ncol(data) > 1) {
        colnames <- paste0(variable, seq(ncol(data)))
      } else {
        colnames <- variable
      }
      if (is.character(dupl)) dupl <- read_binary(paste0(folder, "/", dupl, ".dat"))
      if (length(dupl) == 1) dupl <- rep(dupl, nrow(data))
      data <- data[mrep(seq(nrow(data)), n = dupl), ] %>% data.frame
      data <- data %>% rename_str(colnames)

    })

  if (!is.null(parnames)) {

    pars <- read_parameters(folder, parnames, combine = combine, flatten = TRUE, as_numeric = as_numeric)
    pars <- pars %>% map_dfc(rep, nrow(data))
    data <- list(data, pars) %>% bind_cols

  }

  return (data)

}
