# Variance decomposition
#=========================
# Here we test our variance decomposition procedure
# We simulate genetic values for a focal locus across a population of individuals
# Genotypes and ecotypes are assigned at random
# We calculate the various variance components within and across ecotypes
# Using our equations and comparing them with brute-force variance calculation

rm(list = ls())

library(tidyverse)

# Simulate the data
gmeans <- c(0.3, 0.5, 0.7)
data <- map_dfc(gmeans, ~ rnorm(100, .x, 0))
colnames(data) <- c("aa", "Aa", "AA")
data <- data %>% pivot_longer(cols = colnames(data), names_to = "genotype")
data <- data %>% mutate(count = str_count(genotype, "A"))
data <- data %>% mutate(ecotype = rbinom(nrow(data), 1, 0.5))

# Summaries per ecotype per genotype
smr <- data %>%
  group_by(ecotype, count) %>%
  summarize(
    n = n(),
    sumgen = sum(value),
    ssqgen = sum(value^2)
  ) %>%
  ungroup()

# Function to extract a matrix of interest
get_matrix <- function(data, variable) {
  data %>% 
    select(ecotype, count, !!!variable) %>%
    pivot_wider(names_from = "count", values_from = variable) %>%
    select(-ecotype) %>%
    as.matrix()
}

# Calculate the marginal sums in a matrix
marginalize <- function(m) {
  m <- cbind(m, apply(m, 1, sum))
  rbind(m, apply(m, 2, sum))
}

# Biased variance estimate
var_raw <- function(x) mean(x^2) - mean(x)^2

# Extract the matrices of interest
gcounts <- marginalize(get_matrix(smr, "n"))
gsumgen <- marginalize(get_matrix(smr, "sumgen"))
gssqgen <- marginalize(get_matrix(smr, "ssqgen"))

# Total population size
n <- gcounts[3, 4]

# Population-wide mean allele count and genetic value  
meang <- gsumgen[3, 4] / n
meanq <- (gcounts[3, 2] + 2 * gcounts[3, 3]) / n 

# Check
meang == mean(data$value)
meanq == mean(data$count)

# Calculate variance in allele counts and covariance between genetic value and
# allele counts
covqg <- (2 * gsumgen[3, 3] + gsumgen[3, 2]) / n - meanq * meang
varq <- (4 * gcounts[3, 3] + gcounts[3, 2]) / n - meanq^2

# Check
covqg == mean(with(data, count * value)) - meanq * meang
varq == var_raw(data$count)

# Calculate the average mutational effect (slope of the regression line)
alpha <- covqg / varq

# Run a linear model
mod <- lm(value ~ count, data = data)
intercept <- mod$coefficients[1]
slope <- mod$coefficients[2]

# Check the slope
alpha == slope

# Compute genotype means and additive expected values
gexpec <- map_dbl(0:2, ~ meang + alpha * (.x - meanq))
gmeans <- map2_dbl(gsumgen[3, 1:3], gcounts[3, 1:3], ~ .x / .y)

# Check
gexpec == map_dbl(0:2, ~ intercept + alpha * .x)
gmeans == tapply(data$value, data$count, mean)

# Compute breeding values and dominance deviations
beta <- gexpec - meang
delta <- gmeans - gexpec

# Compute additive variance and check
varA <- apply(gcounts, 1, function(x) sum(x[1:3] * beta^2) / x[4] - (sum(x[1:3] * beta) / x[4])^2)
rbind(varA, apply(gcounts, 1, function(x) var_raw(do.call("c", map2(beta, x[1:3], ~ rep(.x, .y))))))

# Compute dominance variance and check
varD <- apply(gcounts, 1, function(x) sum(x[1:3] * delta^2) / x[4] - (sum(x[1:3] * delta) / x[4])^2)
rbind(varD, apply(gcounts, 1, function(x) var_raw(do.call("c", map2(delta, x[1:3], ~ rep(.x, .y))))))

# Function to compute a given residual variance component (interaction or all non-additive)
get_resid_var <- function(i, from) {
  
  msq <- gssqgen[i, 4] - 2 * sum(from * gsumgen[i, 1:3]) + sum(gcounts[i, 1:3] * from^2) 
  msq <- msq / gcounts[i, 4]
  sqm <- gsumgen[i, 4] - sum(gcounts[i, 1:3] * from)
  sqm <- (sqm / gcounts[i, 4])^2
  msq - sqm
  
}

# Compute the interaction and non-additive variance
varI <- map_dbl(1:3, get_resid_var, gmeans)
varN <- map_dbl(1:3, get_resid_var, gexpec)

# Add interaction and non-additive residuals to the individual-wise dataset
data <- data %>%
  mutate(
    epsilon = value - gmeans[count + 1],
    gamma = value - gexpec[count + 1]
  )

# Check
rbind(varI, c(tapply(data$epsilon, data$ecotype, var_raw), var_raw(data$epsilon)))
rbind(varN, c(tapply(data$gamma, data$ecotype, var_raw), var_raw(data$gamma)))

# Plot ther regression
ggplot(data, aes(x = count, y = value, color = ecotype)) +
  geom_point() +
  geom_abline(intercept = intercept, slope = alpha)

# Conclusion here: the equations used in variance decomposition are valid
# indeed the non-additive variance is zero when the simulation is entirely additive
# this means there must be a problem in the C++ implementation