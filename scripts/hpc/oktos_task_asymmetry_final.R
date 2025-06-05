# Usage: Rscript oktos_task_asymmetry_finalmodel.R <trial_> <model_index> <working_dir> e.g.
# Usage: Call oktos_task_final_asymmetry_.pbs

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
if (!grepl("^[0-1]$", args[2])) {
  stop("Model index must be a number between 0 and 1")
}
model_index <- as.integer(args[2])

# Define trial type
trial_type_ <- ifelse(task_ == "AO", "AO_target-tones", "AXCPT_AB-target")

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

cat("Fitting asymmetry models for ", trial_type_," across regions.\n")

base_path <- "./" # change if needed

# Load the data
model_data <- read.csv(paste0(base_path, task_, "_asymmetry.csv")) %>%
  filter(condition == trial_type_) %>%
  mutate(
    timepoint = factor(timepoint, levels = c("ses-01", "ses-02", "ses-03"), 
                      labels = c("BAS", "POST", "FUP")),
    group = factor(group, levels = c(0, 1), labels = c("NR", "RESP"))
  )

# 3. Transform applied (0-1) to allow for zero-one inflated beta model.
model_data$asymmetry_transformed <- (model_data$asymmetry + 1) / 2

## ---- Specify model parameters (not required but advised)  -------------------------------------------------
iterations <- 7000 # No of iterations
warmup <- 2000 # Warmup samples to be discarded
init <- 0.5
control <- list(
  adapt_engaged = TRUE,
  adapt_delta = 0.999, # increased from default of 0.8
  stepsize = 0.05, # 0.05 default
  max_treedepth = 15
)

# Priors
priors <- 
    c(
        # Priors for intercept, and population effects
        prior(normal(0, 0.25), class = b),  # Sets for factors and their interaction. Try 0, 1
        prior(normal(1, 0.5), class = b, coef = Intercept), # Sets for intercept specifically. I think this is too low - Try 0, 1.5
        
        # Priors for random effects
        prior(gamma(60, 0.1), class = phi),
        prior(exponential(3), class = sd),
        prior(exponential(3), class = sd, dpar = zoi),  # Shrinks variance
        prior(exponential(3), class = sd, dpar = coi),  # Shrinks variance
        
        # 
        prior(normal(-0.5, 0.5), class = Intercept, dpar = zoi), # Weakly shrinks towards fewer zeros
        prior(normal(-0.5, 0.5), class = Intercept, dpar = coi), # Weakly shrinks towards fewer zeros. Try (-0.5, 1)
        
        prior(normal(0, 0.25), class = b, dpar = zoi),
        prior(normal(0, 0.25), class = b, dpar = coi) 
        
    )

# Define zero-one-inflated-beta model formulas
formula <- bf(
  asymmetry_transformed ~ 0 + Intercept + group * timepoint + (1 + timepoint + group || subject:regions),
  zoi ~ group * timepoint + (1 | subject:regions), 
  coi ~ group * timepoint + (1 | subject:regions),
  family = zero_one_inflated_beta()
)

# Conditional model selection based on model_index
if (model_index == 0) {
  model <- brm(
    formula,
    data = model_data,
    chains = 4,
    cores = 4,
    iter = iterations,
    warmup = warmup,
    prior = priors,
    init = init,
    sample_prior = "only",
    save_pars = save_pars(all = TRUE),
    seed = 22,
    file = paste0(base_path, trial_type_, "_asymmetry_finalmodel_prior"),
    control = control
  )
} else if (model_index == 1) {
  model <- brm(
    formula,
    data = model_data,
    chains = 4,
    cores = 4,
    iter = iterations,
    warmup = warmup,
    prior = priors,
    init = init,
    sample_prior = "yes",
    save_pars = save_pars(all = TRUE),
    seed = 22,
    file = paste0(base_path, trial_type_, "_asymmetry_finalmodel"),
    control = control
  )
} else {
  stop("Invalid model index. Please provide a model index between 0 and 1.")
}