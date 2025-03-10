# Usage: Rscript oktos_brms_local_finalmodel.R <metric> <model_index> <working_dir> e.g.
# PBS: Run oktos_local_finalmodel_hpc.pbs

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
library(tidybayes)
library(tidyverse)
library(readxl)
library(furrr)
library(lubridate)
library(dplyr)

cat("Fitting final models for the", metric_, " across channels.\n")

base_path <- "./" # change if needed

# 1. Load the data
model_data <- read.csv(paste0(base_path, "OKTOS_te_local.csv")) %>%
  filter(measure == metric_) %>%
  mutate(
    timepoint = factor(timepoint, levels = c("ses-01", "ses-02", "ses-03"), labels = c("BAS", "POST", "FUP")),
    group = factor(group, levels = c(0, 1), labels = c("NR", "RESP")),
    condition = factor(condition),
    region = relevel(factor(region), "Central_Midline")
  )

# 2. Rescale the DV
model_data$value.s <- scale(model_data$value)[, 1]

### ANALYSIS SECTION ###
## ---- Specify model parameters (not required but advised)  -------------------------------------------------
iterations <- 7000      # No of iterations
warmup <- 2000 # Warmup samples to be discarded
init <- 0
control <- list(
  adapt_engaged = TRUE,
  adapt_delta = 0.999, #increased from default of 0.8
  stepsize = 0.05, # 0.05 default
  max_treedepth = 15
)

## ---- Set priors -------------------------------------------------
prior <- c(
  # Priors for intercept, and population effects
  prior(normal(0, 0.5), class = "b"), # Sets for factors and their interaction
  prior(student_t(3, 0, 1), class = "b", coef = "Intercept"), # Sets for intercept specifically

  # Priors for random effects
  prior(normal(0, 1), class = "sd"),
  prior(lkj(2), class = "cor"),

  # Prior for residual standard deviation
  prior(exponential(2), class = "sigma")
)

# Conditional model selection based on model_index
if (model_index == 0) {
  model <- brm(
    value.s ~ 0 + Intercept + group * timepoint * condition + (1 + timepoint + condition + group | subject:channels),
    data = model_data,
    family = ifelse(metric_ == "Btwn", "exgaussian", 
                  ifelse(metric_ == "ClCoef", "student", "skew_normal")),
    chains = 4,
    cores = 4,
    iter = iterations,
    warmup = warmup,
    prior = prior,
    init = init,
    sample_prior = "only",
    save_pars = save_pars(all = TRUE),
    seed = 22,
    file = paste0(base_path, metric_, "_finalmodel_prior"),
    control = control
  )

} else if (model_index == 1) {
  model <- brm(
    value.s ~ 0 + Intercept + group * timepoint * condition + (1 + timepoint + condition + group | subject:channels),
    data = model_data,
    family = ifelse(metric_ == "Btwn", "exgaussian", 
                  ifelse(metric_ == "ClCoef", "student", "skew_normal")),
    chains = 4,
    cores = 4,
    iter = iterations,
    warmup = warmup,
    prior = prior,
    init = init,
    sample_prior = "yes",
    save_pars = save_pars(all = TRUE),
    seed = 22,
    file = paste0(base_path, metric_, "_finalmodel"),
    control = control
  )
} else {
  stop("Invalid model index. Please provide a model index between 0 and 1.")
}

cat(paste0("Final model fitting for", metric_, " completed for model index:", model_index, "\n"))