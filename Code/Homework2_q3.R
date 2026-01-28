## ---------------------------
##
## Script name: Homework2_q3.R
##
## Author: Trent VanHawkins
##
## Date Created: 2026-01-24
##
##
## ---------------------------
## Load up libraries we need
library(survival)
library(tidyverse)
library(here)
# Map to proper parameterization ------------------------------------------
## Shape/Rate
alpha <- 2
beta <- 0.5

## Shape/Scale
sigma <- 1/beta^(1/alpha)

## Exponential Rate
lambda <- 1/3

# Simulate Failure and Censoring times ------------------------------------
set.seed(404)
n <- 1000
x <- rweibull(n, alpha, sigma)
c <- rexp(n, lambda)

# Define delta and t ------------------------------------------------------
t <- ifelse(x < c, x, c)
delta <- as.numeric(t == x)


# Fit survival function ---------------------------------------------------
fit <- survreg(Surv(t, delta, type = "right") ~ 1, dist="weibull")

## Extract estimated parameters
est_intercept <- coef(fit)
est_scale_tau <- fit$scale

## Convert them back to your alpha and beta
est_alpha <- 1 / est_scale_tau
est_sigma <- as.numeric(exp(coef(fit)))
est_beta  <- exp(-est_intercept / est_scale_tau)

# Compute the predicted values from some values of t ---------------------------
plot_times <- seq(min(t), max(t), length.out = 100)

## Compute true and predicted values
S_t   <- pweibull(plot_times, shape = alpha, scale = sigma, lower.tail = FALSE)
S_hat <- pweibull(plot_times, shape = est_alpha, scale = est_sigma, lower.tail = FALSE)


# Compute the variance and CI ---------------------------------------------
get_ci <- function(t, fit) {
  mu <- as.numeric(coef(fit))
  tau <- fit$scale
  
  # Calculate Z at time t
  z <- (log(t) - mu) / tau
  
  # Gradient (already on log-log scale)
  grad_mu     <- -1 / tau
  grad_logtau <- -z
  grad_vec    <- c(grad_mu, grad_logtau)
  
  # Variance on the log-log scale
  var <- drop(t(grad_vec) %*% vcov(fit) %*% grad_vec)
  se  <- sqrt(var)
  
  # Confidence Interval
  z_lower <- z - 1.96 * se
  z_upper <- z + 1.96 * se
  
  # Back-transform to Survival (0,1) scale
  s_lower <- exp(-exp(z_upper))
  s_upper <- exp(-exp(z_lower))
  
  return(c(lower = s_lower, upper = s_upper))
}

CI <- sapply(plot_times, function(t) get_ci(t, fit))

## Create a dataframe for plotting
res_df <- data.frame(
  time = plot_times,
  True = S_t,
  MLE = S_hat,
  lower_ci = CI[1,],
  upper_ci = CI[2,]
) %>% 
  pivot_longer(cols = c(True, MLE), names_to = "type", values_to = "survival")

colors <- RColorBrewer::brewer.pal(3, "Set2")
## Plot it
ggplot(res_df, aes(x = time, y = survival, color = type)) +
  geom_ribbon(data = res_df %>% filter(type == "MLE"), 
              aes(ymin = lower_ci, ymax = upper_ci, x = time),
              fill = colors[2], 
              alpha = 0.2, 
              inherit.aes = FALSE) + 
  geom_line(size = 1) +
  theme_minimal() +
  labs(title = "True vs. Estimated Weibull Survival Function",
       subtitle = "95% Confidence Interval (Log-Log Scale)",
       y = "S(t)", 
       x = "Time", 
       color = "Type") +
  scale_color_manual(values = c("MLE" = colors[2], "True" = colors[3]))

ggsave("Weibull_survplot.png", path = here("Homework/Figures/"), width = 6, height = 4, units = "in")
