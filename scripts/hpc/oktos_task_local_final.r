# Usage: Rscript oktos_brms_btwn_modelcomp.R <metric> <model_index> <working_dir> e.g.
# PBS: Run oktos_local_hpc.pbs

# Get the metric name from command-line argument
# Validate arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Please provide task, model index, and working directory.")
}

task_ <- args[1]
if (!task_ %in% c("AO", "AXCPT")) {  # Add valid task validation
  stop("Task must be either 'AO' or 'AX'")
}

metric_ <- args[2]

# Validate model index is numeric
if (!grepl("^[0-1]$", args[3])) {
  stop("Model index must be a number between 0 and 1")
}
model_index <- as.integer(args[3])

# Define trial type
trial_type_ <- ifelse(task_ == "AO", "AO_target-tones", "AXCPT_AB-target")

# Working dir
cat("working dir ", getwd(), "\n")
# check if args[3] is a valid directory
working_dir <- args[4]
if (!dir.exists(working_dir)) {
  stop(paste0("The working directory (", working_dir, ") provided does not exist."))
}
setwd(working_dir)
cat("Working dir set to", getwd(), "\n")

library(brms)
library(tidybayes)
library(tidyverse)
library(readxl)
library(furrr)
library(lubridate)
library(dplyr)

cat("Fitting asymmetry models for ", trial_type_," across regions.\n")

base_path <- "./" # change if needed

# 1. Load the data
model_data <- read.csv(paste0(base_path, task_, "_local.csv")) %>%
  filter(trial_type == trial_type_ & measure == metric_) %>%
  mutate(
    timepoint = factor(timepoint, levels = c("ses-01", "ses-02", "ses-03"), labels = c("BAS", "POST", "FUP")),
    group = factor(group, levels = c(0, 1), labels = c("NR", "RESP")),
    channels = relevel(factor(channels), "Cz")
  )

## ---- Specify model parameters (not required but advised)  -------------------------------------------------
iterations <- 7000 # No of iterations
init <- 0.1
warmup <- iterations / 2 # Warmup samples to be discarded
control <- list(
  adapt_engaged = TRUE,
  adapt_delta = 0.999, # increased from default of 0.8
  stepsize = 0.05, # 0.05 default
  max_treedepth = 15
)

# Specify formula conditionally
formula <- if (metric_ %in% c("Btwn", "ClCoef")) {
  bf(
    value ~ 0 + Intercept + group * timepoint + (1 + timepoint + group|| subject:channels),
    hu ~ group * timepoint + (1|| subject:channels)
  )
} else {
  bf(
    value ~ 0 + Intercept + group * timepoint + (1 + timepoint + group || subject:channels),  # returns coefficients for the beta distribution’s mean parameter (mu)
    zi ~ group * timepoint + (1|| subject:channels)
    )
}

priors <- if (metric_ %in% c("Btwn", "ClCoef")) {
    c(
    # Priors for intercept, and population effects
    prior(normal(0, 0.25), class = b),  # Sets for factors and their interaction
    prior(normal(-2, 0.5), class = b, coef = Intercept), # Sets for intercept specifically

    # Priors for random effects
    prior(exponential(3), class = sd),
    prior(exponential(3), class = sd, dpar = hu),  # Shrinks variance

    prior(normal(-0.5, 0.5), class = Intercept, dpar = hu), # Weakly shrinks towards fewer zeros
    prior(normal(0, 0.25), class = b, dpar = hu) # Weakly shrinks towards fewer zeros,
    )
} else {
    c(
    # Priors for intercept, and population effects
    prior(normal(0, 0.25), class = b),  # Sets for factors and their interaction
    prior(normal(-2, 0.5), class = b, coef = Intercept), # Sets for intercept specifically

    # Priors for random effects
    prior(exponential(3), class = sd),
    prior(exponential(3), class = sd, dpar = zi),  # Shrinks variance

    prior(normal(-0.5, 0.5), class = Intercept, dpar = zi), # Weakly shrinks towards fewer zeros
    prior(normal(0, 0.25), class = b, dpar = zi) # Weakly shrinks towards fewer zeros,
    )
}

# Specify the family conditionally
model_family <- ifelse(metric_ %in% c("Btwn", "ClCoef"), "hurdle_gamma", "zero_inflated_negbinomial")

# Conditional execution based on model_index
if (model_index == 0) {
  model <- brm(
    formula,
    data = model_data,
    family = model_family,
    chains = 4,
    cores = 4,
    prior = priors,
    init = init, 
    iter = iterations,
    warmup = warmup,
    sample_prior = "only",
    save_pars = save_pars(all = TRUE),
    seed = 22,
    file = paste0(base_path, task_, "_", trial_type_, "_", metric_, "_finalmodel_prior"),
    control = control
  )
} else if (model_index == 1) {
  model <- brm(
    formula,
    data = model_data,
    family = model_family,
    chains = 4,
    cores = 4,
    prior = priors,
    init = init, 
    iter = iterations,
    warmup = warmup,
    sample_prior = "yes",
    save_pars = save_pars(all = TRUE),
    seed = 22,
    file = paste0(base_path, task_, "_", trial_type_, "_", metric_, "_finalmodel"),
    control = control
  )
} else {
  stop("Invalid model index. Please provide a model index between 0 and 1.")
}

cat("Final model fitting completed for model index:", model_index, "\n")