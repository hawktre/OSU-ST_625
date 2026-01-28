library(ggplot2)
library(dplyr)
library(here)

hazard <- function(t, alpha, beta){
  alpha/(beta + t)
}

t <- seq(1, 20, 0.1)
alpha = seq(1, 3, 1)
beta = seq(1, 3, 1)

hazard_grid = expand.grid(t = t, alpha = alpha, beta = beta)

hazard_grid$hazard <- apply(hazard_grid, 1, function(x) {
  hazard(x["t"], x["alpha"], x["beta"])
})

hazard_grid |> 
  mutate(param_string = factor(paste0("alpha = ", alpha, ", beta = ", beta))) |> 
  ggplot(aes(x = t, y = hazard, color = param_string)) + 
  geom_line()

hazard_grid |> 
  mutate(
    alpha = factor(alpha),
    beta  = factor(beta)
  ) |> 
  ggplot(aes(x = t, y = hazard, color = beta, group = beta)) + 
  geom_line(linewidth = 1) +
  facet_wrap(
    ~ alpha, 
    nrow = 1,
    labeller = labeller(alpha = ~ paste0("alpha = ", .x))
  ) +
  labs(
    x = "t",
    y = expression(lambda(t) == alpha / (beta + t)),
    color = expression(beta),
  ) +
  theme_bw() + 
  theme(legend.position = "bottom")

ggsave("hazard_plot.png", path = here("Homework/Figures"), units = "in", width = 8, height = 4)