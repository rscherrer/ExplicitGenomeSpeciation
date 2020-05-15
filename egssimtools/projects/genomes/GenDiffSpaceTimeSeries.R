rm(list = ls())

library(egssimtools)
library(tidyverse)
library(ggsim)
library(patchwork)

# Read the data
root <- "/media/raphael/bigass/simulations/EGS/genomes/"
data <- readRDS(paste0(root, "simulations.rds"))
backup <- data
data <- backup

head(data)

# Function to plot one pairwise scatterplot
plot_this <- function(data, x, y, hsymmetry = TRUE, shape = 21) {

  final <- data %>% filter(time == max(data$time))

  p <- data %>%
    ggplot(data = data, mapping = aes(x = get(x), y = get(y))) +
    aes(color = scaleI1, alpha = seed, fill = scaleI1) +
    geom_path()
  if (hsymmetry) {
    p <- p + geom_point(data = final, color = "black", size = 4, alpha = 1)
  } else {
    p <- p + geom_point(
      data = final, color = "black", size = 4, alpha = 1, shape = shape
    )
  }
  p <- p +
    scale_color_manual(
      name = expression(sigma[I]), values = c("firebrick1", "deepskyblue4"),
      labels =  c("0", "1")
    ) +
    scale_fill_manual(
      name = expression(sigma[I]), values = c("firebrick1", "deepskyblue4"),
      labels =  c("0", "1")
    ) +
    scale_alpha_manual(
      values = runif(nlevels(factor(data$seed)), min = 0.79, max = 0.81)
    ) +
    theme_bw() +
    guides(alpha = FALSE, fill = FALSE) +
    labs(x = x, y = y)

  if (hsymmetry) {
    p <- p +
      aes(shape = hsymmetry, lty = hsymmetry) +
      scale_shape_manual(name = "h", values = c(21, 24), labels = c("0", "1")) +
      scale_linetype_manual(name = "h", values = c(1, 2), labels = c("0", "1"))
  }

  p

}

# Function to produce a genetic differentiation space figure (for a given architecture)
combine_plots <- function(data, hsymmetry = TRUE, shape = 21) {

  # Rearrange the data to wide-format
  cols <- c("seed", "scaleI1", "location", "time", "genome_Fst")
  if (hsymmetry) cols <- c(cols, "hsymmetry")
  data <- data %>%
    select_at(cols) %>%
    spread_(key_col = "location", value_col = "genome_Fst")

  # Perform PCA
  pca <- data[, (length(cols) - 1):ncol(data)] %>%
    as.matrix %>%
    prcomp(center = TRUE, scale = TRUE)

  # Check out the variance explained
  scree <- tibble(sdev = pca$sdev) %>%
    mutate(pc = 1:nrow(.)) %>%
    ggplot(aes(x = pc, y = sdev)) +
    geom_bar(stat = "identity") +
    theme_bw() +
    labs(x = "Principal component", y = "Standard deviation") +
    ggtitle("Scree plot")
  scree

  # The components are very even, suggesting this architecture
  # varies in Fst along many different combinations of genes
  # But do simulations cluster in this space?

  data <- cbind(data, pca$x %>% data.frame)

  # Make each pairwise scatterplot in a row
  pairs <- 1:5
  pairs <- expand.grid(pairs, pairs)
  pairs <- pairs[pairs[, 1] != pairs[, 2], ]
  pairs[pairs[, 1] > pairs[, 2], ] <- pairs[pairs[, 1] > pairs[, 2], c(2, 1)]
  pairs[, 1] <- paste0("PC", pairs[, 1])
  pairs[, 2] <- paste0("PC", pairs[, 2])
  pairs <- pairs %>% distinct()
  p <- map2(pairs[, 1], pairs[, 2], ~ plot_this(data, .x, .y, hsymmetry, shape))

  # Assemble the plots with patchwork
  fig <- p %>% wrap_plots() +
    plot_layout(guides = "collect") +
    plot_annotation(title = "Genetic differentiation space")

  fig <- fig / scree + plot_layout(heights = c(4, 1))
  fig

}

# Split by architecture and apply the plotting function
newdata <- data %>%
  group_by(archfile) %>%
  nest() %>%
  mutate(figure = map(data, combine_plots))

# Check one plot
newdata$figure[[1]]

# Saving function
save_this <- function(x, y) {
  figname <- sprintf("gendiffspace_timeseries_%s.png", gsub(".txt", "", x))
  ggsave(figname, y, width = 8, height = 7, dpi = 400)
}

# Save the figures
with(newdata, walk2(archfile, figure, save_this))



# Further split by ecological mode
newdata2 <- data %>%
  group_by(archfile, hsymmetry) %>%
  nest() %>%
  mutate(shape = ifelse(hsymmetry == 0, 21, 24)) %>%
  mutate(figure = map2(data, shape, ~ combine_plots(.x, hsymmetry = FALSE, shape = .y)))

newdata2
newdata2$figure[[2]]

figname <- "gendiffspace_timeseries_%s_h_%s.png"
newdata2 <- newdata2 %>%
  mutate(figname = sprintf(figname, gsub(".txt", "", archfile), hsymmetry))

# Saving function
save_this2 <- function(figure, figname) {
  ggsave(figname, figure, width = 8, height = 7, dpi = 400)
}

# Save the figures
with(newdata2, walk2(figure, figname, save_this2))
