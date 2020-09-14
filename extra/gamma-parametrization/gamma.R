rm(list = ls())

library(tidyverse)

# Bilateral Gamma distribution
rbigamma <- function(n, shape = 1, scale = 1) {

  x <- rgamma(n, shape, scale = scale)
  to_flip <- as.logical(rbinom(n, size = 1, prob = 0.5))
  x[to_flip] <- -x[to_flip]
  return(x)

}

# Sample additive genetic effects
sample_add <- function(n, shape = 1, scale = 1) {

  # Sample independent locus effect sizes from a bilateral Gamma
  eta <- rbigamma(n, shape, scale)
  eta <- eta / sum(eta^2)

  # Sample independent gene expression levels (= genotypes)
  xi <- sample(-1:1, size = n, replace = TRUE)

  # Get distribution of independent effects
  return(eta * xi)

}

# Sample epistatic interaction effects
sample_epi <- function(n, shape = 1, scale = 1) {

  # Sample interaction weights from a bilateral Gamma
  omega <- rbigamma(n, shape, scale)
  omega <- omega / sum(omega^2)

  # Sample pairs of gene expression levels
  xi_i <- sample(-1:1, size = n, replace = TRUE)
  xi_j <- sample(-1:1, size = n, replace = TRUE)

  # Get distribution of interaction effects across the genome
  return(omega * xi_i * xi_j)

}

phi <- sample_add(10000, 5, 1)
psi <- sample_epi(100000, 5, 1)

hist(phi)
hist(psi)

n <- 1000
e <- 10000
shape_add <- 1
shape_epi <- 0.12
scale_add <- 1
scale_epi <- 1

set.seed(42)
x <- purrr::map_dbl(seq(10000), ~ sum(sample_add(n, shape_add, scale_add)))
y <- purrr::map_dbl(seq(10000), ~ sum(sample_epi(e, shape_epi, scale_epi)))

data <- tibble(x = x, y = y)
data <- data %>%
  pivot_longer(cols = c("x", "y"), names_to = "variable", values_to = "value")

ggplot(data, aes(x = value, fill = variable)) +
  geom_density(alpha = 0.3) +
  theme_bw() +
  scale_fill_manual(
    labels = c("Additive", "Epistatic"),
    values = c("deepskyblue4", "firebrick1")
  ) +
  labs(x = "Trait value", y = "Density", fill = NULL)

#ggsave("gamma_param.png", width = 5, height = 3, dpi = 300)
