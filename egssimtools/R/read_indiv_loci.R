#' Read individual loci
#'
#' Makes a data frame with alleles at each locus for each individual.
#' Can be very large
#'
#' @param root Path to the simulation folder
#' @param nloci Number of loci. Guessed using `guess_nloci` if NULL.
#' @param t Optional time points to read. All if NULL.
#' @param extra Optional names of extra individual-wise variables to attach to
#' the output data frame (read using `read_pop`)
#' @param freezerfile File name of the freezer file
#' @param locifile File name of the loci values file
#' @param architecture Whether to add components of the genetic architecture
#' @param archfile Genetic architecture file name
#'
#' @return A data frame
#'
#' @examples
#'
#' \dontrun{
#'
#' root <- "data/example_1"
#'
#' read_indiv_loci(root)
#'
#' }
#'
#' @export

read_indiv_loci <- function(
  root,
  nloci = NULL,
  t = NULL,
  extra = NULL,
  freezerfile = "freezer.dat",
  locifile = "locivalues.dat",
  architecture = TRUE,
  archfile = "architecture.txt"
) {

  filename <- file.path(root, freezerfile)
  genome <- read_bitset(filename)

  if (is.null(nloci)) nloci <- guess_nloci(root)

  # Number of zeros added at the end of each individual genome
  nextra <- 64 - ((2 * nloci) %% 64)

  # Number of bits per individual genome
  nbits <- (2 * nloci) + nextra

  # Number of individuals
  ninds <- length(genome) / nbits

  # Make a large dataset
  inds <- rep(seq(ninds), each = nbits)
  loci <- rep(c(rep(seq(nloci), 2), rep(NA, nextra)), ninds)
  haps <- rep(c(rep(seq(2), each = nloci), rep(NA, nextra)), ninds)
  data <- tibble::tibble(
    allele = genome,
    locus = loci,
    ind = inds,
    hap = haps
  )

  # Remove the extra zeros (they were added to convert individual genomes
  # into 64bit integers)
  data <- data %>% tidyr::drop_na()

  # Convert to the wide format (don't thinks we will need the long-format
  # in this study) = locus-wise
  data <- data %>% dplyr::mutate(hap = stringr::str_replace(hap, "^", "hap"))
  data <- data %>% tidyr::pivot_wider(names_from = "hap", values_from = "allele")

  # Read the time point when each individual lived
  popsizes <- read_sim(root, "population_size")
  if (is.null(t)) t <- unique(popsizes$time)
  popsizes <- popsizes %>% dplyr::filter(time %in% t)
  times <- do.call("c", with(popsizes, purrr::map2(time, population_size, ~ rep(.x, .y))))
  data$time <- rep(times, each = nloci)

  # Derive allele counts from the alleles
  data <- data %>% dplyr::mutate(allcount = hap1 + hap2)

  # Add loci values
  data$value <- read_binary(file.path(root, locifile))

  # Add additional individual data if needed
  if (!is.null(extra)) {
    extradata <- read_pop(root, extra)
    extradata <- extradata %>% dplyr::filter(time %in% t)
    extradata <- extradata %>% dplyr::rename_all(~ stringr::str_replace(.x, "individual_", ""))
    extradata <- extradata %>% dplyr::mutate(ind = seq(nrow(extradata)))
    data <- data %>% dplyr::right_join(extradata)
  }

  # Add genetic architecture components if needed
  if (architecture) {
    arch <- read_arch_genome(root, archfile)
    data <- data %>% dplyr::right_join(arch)
  }

  return(data)

}
