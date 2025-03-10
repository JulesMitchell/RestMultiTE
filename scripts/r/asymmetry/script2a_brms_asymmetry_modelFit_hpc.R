# Usage: Rscript oktos_brms_asymmetry_finalmodel.R <metric> <working_dir> e.g.
# Usage: Call asymmetry_final_hpc.pbs

# Get the metric name from command-line argument
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
    stop("Please provide model index, metric, and working directory.")
}

metric_ <- args[1]
model_index <- as.integer(args[2])  # The model index passed from PBS script (0 to 5)
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

#3. No scaling applied as -1 to 1 is meaningful.

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

# Priors
prior = c(
  # Priors for intercept, and population effects
  prior(normal(0, 0.1), class = 'b'),  # Sets for factors and their interaction
  prior(normal(0, 0.25), class = 'b', coef='Intercept'), # Sets for intercept specifically

  # Priors for random effects
  prior(normal(0, 0.25), class = 'sd'),
  prior(lkj(2), class='cor'),

  # Prior for residual standard deviation
  prior(exponential(2), class = "sigma")
)

# Conditional model selection based on model_index
if (model_index == 0) {
  finalModel_prior <- brm(asymmetry| resp_trunc(lb = -1, ub=1) ~ 0 + Intercept + group * timepoint * condition + (1 + timepoint + condition| subject:regions),
            data = model_data, 
            family = student(), 
            chains = 4,  # Number of chains
            cores = 4,   # Number of CPU cores
            iter = iterations,  # Number of iterations
            warmup = warmup,  # Warmup iterations
            prior = prior,
            init = init,
            sample_prior = 'only',
            save_pars = save_pars(all = TRUE),
            seed = 22,
            file = paste0(base_path, metric_, "_finalmodel_prior"), # Uncomment save model object
            control = control # adjust if divergence issues emerge
)
  
} else if (model_index == 1) {
  finalModel <- brm(asymmetry| resp_trunc(lb = -1, ub=1) ~ 0 + Intercept + group * timepoint * condition + (1 + timepoint + condition| subject:regions),
            data = model_data, 
            family = student(),  # Gaussian likelihood
            chains = 4,  # Number of chains
            cores = 4,   # Number of CPU cores
            iter = iterations,  # Number of iterations
            warmup = warmup,  # Warmup iterations
            prior = prior,
            init = init,
            sample_prior = 'yes',
            save_pars = save_pars(all = TRUE),
            seed = 22,
            file = paste0(base_path, metric_, "_finalmodel"),
            control = control # adjust if divergence issues emerge
)

} else {
  stop("Invalid model index. Please provide a model index between 0 and 1.")
}

cat("Final model asymmetry fitting completed for model index:", model_index, "\n")