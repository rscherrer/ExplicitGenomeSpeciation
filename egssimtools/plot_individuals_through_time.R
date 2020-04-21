rm(list = ls())

simulation <- "../build/release"

# Plot a simulation through time
data <- read_individuals(simulation, "individual_trait")
data <- lapply(data, function(data) sapply(data, "[", 1))
data <- data.frame(time = mrep(as.numeric(names(data)), sapply(data, length)), data = do.call("c", data))

p <- ggplot(data, aes(x = time, y = data)) +
  geom_bin2d(bins = 15) +
  theme_bw() +
  scale_fill_continuous(type = "viridis") +
  xlab("Time (generations)") +
  ylab("Ecological trait") +
  labs(fill = "Count")
p
