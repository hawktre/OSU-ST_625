## ---------------------------
##
## Script name: Homework4_Q1.R
##
## Purpose of script:
##
## Author: Trent VanHawkins
##
## Date Created: 2026-02-13
## ---------------------------

## view outputs in non-scientific notation

options(scipen = 6, digits = 4) 

## ---------------------------

## load up the packages we will need:  (uncomment as required)

library(tidyverse)
library(here)
library(KMsurv)
library(survival)

data(kidney)
head(kidney)


# Fit the model -----------------------------------------------------------

perc_patients <- kidney[kidney$type == 2,]
fit <- survreg(Surv(time, delta) ~ 1, data = perc_patients )
summary(fit)

## Extract estimated parameters
est_intercept <- coef(fit)
est_scale_tau <- fit$scale
