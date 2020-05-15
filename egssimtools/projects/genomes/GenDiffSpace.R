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

# Downsample to the last time point
data <- data %>% filter(time == 19900)

# Function for a given pairwise scatterplot
plot_this <- function(data, x, y) {

  #library(plyr)
  #find_hull <- function(df, x, y) df[chull(df[[x]], df[[y]]),]
  #data <- data %>% mutate(group = paste(hsymmetry, scaleI1))
  #hulls <- ddply(data, "group", find_hull, x, y)

  data %>%
    ggplot(data = data, mapping = aes(x = get(x), y = get(y))) +
    aes(fill = scaleI1, shape = hsymmetry, color = scaleI1, alpha = hsymmetry) +
    geom_point(size = 3, color = "black") +
    stat_ellipse() +
    #geom_polygon(data = hulls) +
    scale_shape_manual(name = "h", values = c(21, 24), labels = c("0", "1")) +
    scale_alpha_manual(name = "h", values = c(0.3, 0.7), labels = c("0", "1")) +
    scale_fill_manual(name = expression(sigma[I]), values = c("firebrick1", "deepskyblue4")) +
    scale_color_manual(name = expression(sigma[I]), values = c("firebrick1", "deepskyblue4")) +
    theme_bw() +
    guides(fill = FALSE) +
    labs(x = x, y = y)

}

# Function to produce a genetic differentiation space figure (for a given architecture)
combine_plots <- function(data) {

  # Rearrange the data to wide-format
  data <- data %>%
    select(seed, hsymmetry, scaleI1, location, genome_Fst) %>%
    spread_(key_col = "location", value_col = "genome_Fst")

  # Perform PCA
  pca <- data[, 4:ncol(data)] %>%
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
  p <- map2(pairs[, 1], pairs[, 2], ~ plot_this(data, .x, .y))

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

# Check plot
newdata$figure[[1]]

# Saving function
save_this <- function(x, y) {
  figname <- sprintf("gendiffspace_%s.png", gsub(".txt", "", x))
  ggsave(figname, y, width = 8, height = 7, dpi = 400)
}

# Save the figures
with(newdata, walk2(archfile, figure, save_this))
