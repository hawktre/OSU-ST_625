## ---------------------------
##
## Script name: Homework3_q1.R
##
## Purpose of script: Create functions to plot the kaplan meier estimator
##
## Author: Trent VanHawkins
##
## Date Created: 2026-01-28
##
##
## ---------------------------
## load up the packages we will need:  (uncomment as required)

library(survival)
library(tidyverse)
library(ggpubr)

aml <- survival::aml

# 1. Hazard Rate Fucntion -------------------------------------------------
## funciton to return lambda_hat values
get_lambda_hat <- function(t, c){
  #get count of failures at each t
  d <- table(t[c == 1])
  failure_times <- as.numeric(names(d)) #times of failures
  y_bar <- sapply(failure_times, function(time) sum(t >= time)) # Number remianing at time t
  
  d/y_bar
}

# Run the function
lambda_hat <- get_lambda_hat(aml$time, aml$status)

# 2. Get the kaplean-meier estimator and plot -----------------------------
failure_times <- as.numeric(names(lambda_hat))

# Get S_km 
S_km <- sapply(failure_times, function(t) prod(1 - lambda_hat[which(failure_times <= t)]))

#Create a data frame for plotting (add a point for time 0)
estimate_df <- data.frame(t = c(0, failure_times, max(aml$time)), S = c(1,S_km, S_km[length(S_km)]))

#Plot it
estimate_df %>% 
  ggplot(aes(x = t, y = S))+
  geom_step()


# 3. log-log confidence bound ---------------------------------------------

loglog_ci <- function(s_km, t, c){
  #Get what we need to construct sigma_hat
  d <- table(t[c == 1]) #number of failures at time t
  failure_times <- as.numeric(names(d)) #tiems of failures
  y_bar <- sapply(failure_times, function(time) sum(t >= time)) # Number remianing at time t
  
  #compute sigma_hat
  inner_sum <- d/(y_bar * (y_bar - d))
  sigma_squared <- sapply(failure_times, function(time) sum(inner_sum[which(failure_times <= time)]))
  sigma_hat <- sqrt(sigma_squared)
  
  ## compute SE on log-log scale
  loglog_se <- 1/abs(log(s_km)) * sigma_hat
  
  ## Compute the confidence interval
  lower_ci <- s_km^(exp(loglog_se * 1.96))
  upper_ci <- s_km^(exp(-loglog_se * 1.96))
  
  return(list(lower = lower_ci, upper = upper_ci))
}

s_ci <- loglog_ci(S_km, aml$time, aml$status)


# 4-5 Show that the plots are the same ------------------------------------
## include confidence bands in data frame
estimate_df$lower_ci <- c(1, s_ci$lower, s_ci$lower[which.min(s_ci$lower)])
estimate_df$upper_ci <- c(1, s_ci$upper, s_ci$upper[which.min(s_ci$upper)])

#Actual sruvival fit and data frame
fit <- survfit(Surv(time, status) ~ 1, data = aml, conf.type = "log-log")
fit_summary <- summary(fit)
survival_df <- data.frame(t = c(0,fit_summary$time,max(aml$time)),
                          S = c(1, fit_summary$surv, min(fit_summary$surv)),
                          lower_ci = c(1, fit_summary$lower, min(fit_summary$lower)),
                          upper_ci = c(1, fit_summary$upper, min(fit_summary$upper)))

colors <- RColorBrewer::brewer.pal(n = 3, name = "Set2")
## Plot both on the same plot
by_hand <- estimate_df %>% 
  ggplot(aes(x = t, y = S))+
  geom_step(color = colors[2], linewidth = 1.25)+
  geom_step(aes(y = lower_ci), linetype = "dashed")+
  geom_step(aes(y = upper_ci), linetype = "dashed")+
  theme_minimal()+
  labs(y = expression(S^KM * (t)),
       x = "Time", 
       title = "Kaplan Meir Estimator (By Hand)")


survival_plot <- survival_df %>% 
  ggplot(aes(x = t, y = S))+
  geom_step(color = colors[3], linewidth = 1.25)+
  geom_step(aes(y = lower_ci), linetype = "dashed")+
  geom_step(aes(y = upper_ci), linetype = "dashed")+
  theme_minimal()+
  labs(y = expression(S^KM * (t)),
       x = "Time", 
       title = "Kaplan Meir Estimator (survfit)")

ggarrange(by_hand, survival_plot, nrow = 1)


# Create log-rank statistic ------------------------------------------------------
## Using survival
lr_survival <- survdiff(Surv(time, status) ~ x, data = aml)

lr_test <- function(t, c, j, group){

  ## Indicators
  failed <- c == 1
  group_ind <- j == group
  ## number of failures
  all_failure_times <- sort(unique(t[failed]))
  d_i <- table(factor(t[failed], levels = all_failure_times)) #number of failures at time t
  d_j <- table(factor(t[failed & group_ind == 1], levels = all_failure_times)) #number of failures at time t in group j
  
  ##Number remaining at time t
  y_bar_i <- sapply(all_failure_times, function(time) sum(t >= time)) # both groups
  y_bar_j <- sapply(all_failure_times, function(time) sum(t[group_ind] >= time)) #group j
  
  # Compute the statistic
  expected_j <- d_i * (y_bar_j / y_bar_i)
  z_j <- sum(d_j - expected_j, na.rm = TRUE)
  
  #Compute the standard error
  zj_var <- sum(d_i*(y_bar_j/y_bar_i)*(1 - (y_bar_j/y_bar_i))*((y_bar_i - d_i)/(y_bar_i - 1)), na.rm = TRUE)
  
  #Compute the p-value
  test_stat <- (z_j/sqrt(zj_var))**2
  pval <- pchisq(abs(test_stat), df = 1, lower.tail = FALSE)
  return(list(Chisq = test_stat, p = pval))
}

# Run the by-hand test
lr_byhand <- lr_test(aml$time, aml$status, aml$x, "Maintained")

#Check if they are the same (without machine-rounding error: 
(lr_byhand$Chisq - lr_survival$chisq) < 0.001
(lr_byhand$p - lr_survival$pvalue) < 0.001
