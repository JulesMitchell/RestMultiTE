# Usage: Rscript oktos_task_sampen_modelcomp.R <task> <model_index> <working_dir> e.g.
# PBS: Run oktos_task_modelcomp_sampen.pbs

# Get the metric name from command-line argument
# Validate arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Please provide task, model index, and working directory.")
}

task_ <- args[1]
if (!task_ %in% c("AO", "AXCPT")) {  # Add valid task validation
  stop("Task must be either 'AO' or 'AX'")
}

# Validate model index is numeric
if (!grepl("^[0-3]$", args[2])) {
  stop("Model index must be a number between 0 and 3")
}
model_index <- as.integer(args[2])

# Define trial type
trial_type_ <- ifelse(task_ == "AO", "target_tones", "AB-target")

# Working dir
cat("working dir ", getwd(), "\n")
# check if args[3] is a valid directory
working_dir <- args[3]
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

cat("Fitting sample entropy models for the", task_, " across channels.\n")

base_path <- "./" # change if needed

# Load and process data
model_data <- read.csv(paste0(base_path, task_, "_entropy_channels.csv")) %>%
  filter(trial_type == trial_type_) %>%
  mutate(
    timepoint = factor(timepoint, levels = c("ses-01", "ses-02", "ses-03"), 
                      labels = c("BAS", "POST", "FUP")),
    group = factor(group, levels = c(0, 1), labels = c("NR", "RESP")),  # Fixed parenthesis
    channel = relevel(factor(channel), "Cz")
  )

## ---- Specify model parameters (not required but advised)  -------------------------------------------------
iterations <- 2000 # No of iterations
warmup <- iterations / 2 # Warmup samples to be discarded
init <- 0.5
control <- list(
  adapt_engaged = TRUE,
  adapt_delta = 0.999, # increased from default of 0.8
  stepsize = 0.05, # 0.05 default
  max_treedepth = 15
)

# Priors
priors <- c(
  # Priors for intercept, and population effects
  prior(normal(0, 0.25), class = b),  # Sets for factors and their interaction
  prior(normal(0, 1), class = b, coef = Intercept), # Sets for intercept specifically

  # Priors for random effects
  prior(exponential(3), class = sd)
)

# Conditional execution based on model_index
if (model_index == 0) {
  model <- brm(
    sampen ~ 0 + Intercept + group * timepoint + (1|| subject:channel),
    data = model_data,
    family = ifelse(task_ == "AO", "Beta", "Beta"),
    chains = 2,
    cores = 4,
    iter = iterations,
    warmup = warmup,
    prior = priors,
    init = init,
    sample_prior = "yes",
    save_pars = save_pars(all = TRUE),
    seed = 22,
    file = paste0(base_path, task_, "_model1"),
    control = control
  )
} else if (model_index == 1) {
  model <- brm(
    sampen ~ 0 + Intercept + group * timepoint + (1 + timepoint|| subject:channel),
    data = model_data,
    family = ifelse(task_ == "AO", "Beta", "Beta"),
    chains = 2,
    cores = 4,
    iter = iterations,
    warmup = warmup,
    prior = priors,
    init = init,
    sample_prior = "yes",
    save_pars = save_pars(all = TRUE),
    seed = 22,
    file = paste0(base_path, task_,  "_model2"),
    control = control
  )
} else if (model_index == 2) {
  model <- brm(
    sampen ~ 0 + Intercept + group * timepoint + (1 + group || subject:channel),
    data = model_data,
    family = ifelse(task_ == "AO", "Beta", "Beta"),
    chains = 2,
    cores = 4,
    iter = iterations,
    warmup = warmup,
    prior = priors,
    init = init,
    sample_prior = "yes",
    save_pars = save_pars(all = TRUE),
    seed = 22,
    file = paste0(base_path, task_,  "_model3"),
    control = control
  )
} else if (model_index == 3) {
  model <- brm(
    sampen ~ 0 + Intercept + group * timepoint + (1 + timepoint + group|| subject:channel),
    data = model_data,
    family = ifelse(task_ == "AO", "Beta", "Beta"),
    chains = 2,
    cores = 4,
    iter = iterations,
    warmup = warmup,
    prior = priors,
    init = init,
    sample_prior = "yes",
    save_pars = save_pars(all = TRUE),
    seed = 22,
    file = paste0(base_path, task_,  "_model4"),
    control = control
  )
} else {
  stop("Invalid model index. Please provide a model index between 0 and 3.")
}

cat("Model fitting completed for model index:", model_index, "\n")