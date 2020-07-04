rm(list = ls())

library(tidyverse)

gmeans <- c(0.2, 0.7, 0.8)
n <- 10
q <- rep(0:2, each = n)
g <- reduce(map(gmeans, ~ rnorm(10, .x, 0.05)), c)
df <- tibble(q = q, g = g, q2 = q + 1)
df_means <- df %>% group_by(q) %>% summarize(g = mean(g))
df_means <- df_means %>% mutate(q2 = q + 1)
df_means <- df_means %>% mutate(xmin = q2 - 0.25, xmax = q2 + 0.25)
df_mean2 <- tibble(q = mean(q), g = mean(g), q2 = q + 1)

mod <- lm(g ~ q2, data = df)
intpt <- mod$coefficients[1]
alpha <- mod$coefficients[2]

df_means$pred <- intpt + alpha * (0:2 + 1)

labels <- c(
  "bar(q)[i]", "bar(g)[i]", "g[ij]^'*'", "tilde(g)[ij]", "g[ijk]", "alpha",
  "beta[ij]", "delta[ij]", "epsilon[ijk]"
)
label <- "a[i]"

p <- ggplot(df, aes(x = factor(q), y = g)) +
  geom_jitter(width = 0.2) +
  geom_abline(slope = alpha, intercept = intpt, linetype = 2) +
  geom_point(data = df_means, mapping = aes(x = q2, y = pred), fill = "gray40", size = 3, shape = 21) +
  geom_point(
    data = df_mean2, mapping = aes(x = q2, y = g),
    size = 5, shape = 21, fill = "darkgrey"
  ) +
  geom_segment(
    data = df_means,
    mapping = aes(x = xmin, xend = xmax, y = g, yend = g), size = 1
  ) +
  annotate("text", x = 2.2, y = 0.5, label = "'('*bar(q)[i]*','~bar(g[i])*')'", parse = TRUE) +
  annotate("text", x = 0.65, y = 0.4, label = "beta[ij]", parse = TRUE) +
  annotate("text", x = 0.65, y = 0.23, label = "delta[ij]", parse = TRUE) +
  annotate("text", x = 0.65, y = 0.06, label = "epsilon[ijk]", parse = TRUE) +
  annotate("text", x = 0.49, y = 0.18, label = "gamma[ijk]", parse = TRUE) +
  annotate("text", x = 1.3, y = min(g), label = "g[ijk]", parse = TRUE, hjust = 0) +
  annotate("text", x = 1.3, y = min(df_means$g), label = "tilde(g)[ij]", parse = TRUE, hjust = 0) +
  annotate("text", x = 1.3, y = min(df_means$pred), label = "g[ij]*'*'", parse = TRUE, hjust = 0) +
  annotate("text", x = 1.8, y = 0.4, label = "'slope'==alpha[i]", parse = TRUE) +
  theme_classic() +
  xlab(parse(text = "'Allele count,'~q[i]")) +
  ylab(parse(text = "'Genetic value,'~g[i]"))
p

ggsave("variance-decomp.pdf", p, width = 4, height = 3)
