# Usage: Rscript oktos_brms_asymmetry_modelcomp.R <metric> <model_index> <working_dir> e.g.

# Get the metric name from command-line argument
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Please provide model index, metric, and working directory.")
}

metric_ <- args[1]
model_index <- as.integer(args[2]) # The model index passed from PBS script (0 to 5)
working_dir <- args[3]

cat("working dir ", getwd(), "\n")
# check if args[3] is a valid directory
if (!dir.exists(working_dir)) {
  stop(paste0("The working directory (", working_dir, ") provided does not exist."))
}
setwd(working_dir)
cat("Working dir set to", getwd(), "\n")

library(brms)
library(ggplot2)
library(tidybayes)
library(tidyverse)
library(readxl)
library(furrr)
library(lubridate)
library(dplyr)

cat("Fitting models for the", metric_, " across regions.\n")

base_path <- "./" # change if needed

# Load the data
model_data <- read.csv(paste0(base_path, "OKTOS_te_asymmetry.csv")) %>%
  mutate(
    timepoint = factor(timepoint, levels = c("ses-01", "ses-02", "ses-03"), labels = c("BAS", "POST", "FUP")),
    group = factor(group, levels = c(0, 1), labels = c("NR", "RESP")),
    condition = factor(condition)
  )

# 3. No scaling applied as -1 to 1 is meaningful.

## ---- Specify model parameters (not required but advised)  -------------------------------------------------
iterations <- 4000 # No of iterations
warmup <- iterations / 2 # Warmup samples to be discarded
init <- 0
control <- list(
  adapt_engaged = TRUE,
  adapt_delta = 0.999, # increased from default of 0.8
  stepsize = 0.05, # 0.05 default
  max_treedepth = 15
)

# Priors
prior_cor <- c(
  # Priors for intercept, and population effects
  prior(normal(0, 0.1), class = "b"), # Sets for factors and their interaction
  prior(normal(0, 0.25), class = "b", coef = "Intercept"), # Sets for intercept specifically

  # Priors for random effects
  prior(normal(0, 0.25), class = "sd"),
  prior(lkj(2), class = "cor"),

  # Prior for residual standard deviation
  prior(exponential(2), class = "sigma")
)

# For models without random slope correlations
prior <- c(
  # Priors for intercept, and population effects
  prior(normal(0, 0.1), class = "b"), # Sets for factors and their interaction
  prior(normal(0, 0.25), class = "b", coef = "Intercept"), # Sets for intercept specifically

  # Priors for random effects
  prior(normal(0, 0.25), class = "sd"),

  # Prior for residual standard deviation
  prior(exponential(2), class = "sigma")
)

# Conditional model selection based on model_index
if (model_index == 0) {
  model <- brm(
    asymmetry | resp_trunc(lb = -1, ub = 1) ~ 0 + Intercept + group * timepoint * condition + (1 | subject:regions),
    data = model_data,
    family = student(),
    chains = 2,
    cores = 4,
    iter = iterations,
    warmup = warmup,
    prior = prior,
    init = init,
    sample_prior = "yes",
    save_pars = save_pars(all = TRUE),
    seed = 22,
    file = paste0(base_path, "asymmetry_model1"),
    control = control
  )
} else if (model_index == 1) {
  model <- brm(
    asymmetry | resp_trunc(lb = -1, ub = 1) ~ 0 + Intercept + group * timepoint * condition + (1 + timepoint| subject:regions),
    data = model_data,
    family = student(),
    chains = 2,
    cores = 4,
    iter = iterations,
    warmup = warmup,
    prior = prior_cor,
    init = init,
    sample_prior = "yes",
    save_pars = save_pars(all = TRUE),
    seed = 22,
    file = paste0(base_path, "asymmetry_model2"),
    control = control
  )
} else if (model_index == 2) {
  model <- brm(
    asymmetry | resp_trunc(lb = -1, ub = 1) ~ 0 + Intercept + group * timepoint * condition + (1 + timepoint + condition | subject:regions),
    data = model_data,
    family = student(),
    chains = 2,
    cores = 4,
    iter = iterations,
    warmup = warmup,
    prior = prior_cor,
    init = init,
    sample_prior = "yes",
    save_pars = save_pars(all = TRUE),
    seed = 22,
    file = paste0(base_path, "asymmetry_model3"),
    control = control
  )
} else if (model_index == 3) {
  model <- brm(
    asymmetry | resp_trunc(lb = -1, ub = 1) ~ 0 + Intercept + group * timepoint * condition + (1 + timepoint + group | subject),
    data = model_data,
    family = student(),
    chains = 2,
    cores = 4,
    iter = iterations,
    warmup = warmup,
    prior = prior_cor,
    init = init,
    sample_prior = "yes",
    save_pars = save_pars(all = TRUE),
    seed = 22,
    file = paste0(base_path, "asymmetry_model4"),
    control = control
  )
} else if (model_index == 4) {
  model <- brm(
    asymmetry | resp_trunc(lb = -1, ub = 1) ~ 0 + Intercept + group * timepoint * condition + (1 + timepoint + condition + group| subject),
    data = model_data,
    family = student(),
    chains = 2,
    cores = 4,
    iter = iterations,
    warmup = warmup,
    prior = prior_cor,
    init = init,
    sample_prior = "yes",
    save_pars = save_pars(all = TRUE),
    seed = 22,
    file = paste0(base_path, "asymmetry_model5"),
    control = control
  )
} else {
  stop("Invalid model index. Please provide a model index between 0 and 5.")
}

cat("Model asymmetry fitting completed for model index:", model_index, "\n")