
smr <- data %>%
  filter(scaleI == "1 0 0", mutation == 0.001) %>%
  mutate(ecosel = fct_rev(factor(ecosel)), hsymmetry = factor(hsymmetry)) %>%
  group_by(hsymmetry, ecosel, dispersal, mutation, scaleA, scaleI, simulation) %>%
  mutate(color = last(SI))

ggplot(smr, aes(x = time, y = EI, alpha = simulation, color = color)) %>%
  facettize(
    smr,
    facet_rows = "ecosel",
    facet_cols = "hsymmetry",
    label_facets = TRUE,
    facet_prefixes = c("s", "h")
  ) +
  geom_line() +
  scale_alpha_manual(values = runif(2000, min = 0.49, max = 0.51)) + # hack
  guides(alpha = FALSE) +
  scale_color_gradient(low = "black", high = "coral")
