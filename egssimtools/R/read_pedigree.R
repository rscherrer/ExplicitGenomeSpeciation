#' Read the results of a pedigree experiment
#'
#' The pedigree experiment saves distributions of offspring of random mating
#' pairs
#'
#' @param root Path to the simulation folder
#' @param k Number of offspring per cross
#' @param tped Frequency of the pedigree experiment, in generations
#' @param tend Number of generations in the simulation
#' @param filename Pedigree file name
#'
#' @return A data frame with information, for each cross, about the male and the female, their ecotypes, their trait values and the trait values of their k offspring. There is an extra column for time point in case the experiment was conducted at multiple time points.
#'
#' @export

read_pedigree <- function(root, k, tped, tend, filename = "pedigree.dat") {

  data <- read_binary(file.path(root, filename))

  by <- 4 + 2 * 3 + k * 3

  data <- data %>%
    split(rep(seq(length(.) / by), each = by)) %>%
    do.call("rbind", .) %>%
    data.frame

  colnames <- c(

    "fem", "mal", "ecof", "ecom",
    do.call("c", map(c("f", "m"), ~ paste0(c("x", "y", "z"), .x))),
    do.call("c", map(seq(k), ~ paste0(c("x", "y", "z"), .x)))

  )

  colnames(data) <- colnames

  t <- seq(tped, by = tped, length.out = tend %/% tped)
  n <- nrow(data) / length(t)

  data$time <- rep(t, each = n)

  return(data)

}
