### MODEL SPECIFICATION AND FITTING###
## ----pkgload, include = FALSE-------------------------------------------------
library(brms)
library(dplyr)
library(xlsx)

## ---- Set working directory as required -------------------------------------------------
outputFile <- "c:/Users/j_m289/Pictures/phd/3. Data Analysis/studies/OKTOS/resting_te/analysis/brms/output/geff_finaloutput.xlsx"

## ----Load Datasets, Clean and Set Factor Levels -------------------------------------------------
#1. Load data
data <- read.csv("c:/Users/j_m289/Pictures/phd/3. Data Analysis/studies/OKTOS/resting_te/analysis/spreadsheets/OKTOS_te_global.csv") %>%
    mutate(
        timepoint = factor(timepoint, levels = c("ses-01", "ses-02", "ses-03"), labels = c("BAS", "POST", "FUP")),
        group = factor(group, levels = c(0, 1), labels = c("NR", "RESP")),
        condition = factor(condition),
      )

# 2. Scaling and plot
# Rescale DV 
model_data$g_ef.s <- scale(model_data$g_ef)[, 1]

### ANALYSIS SECTION ###
## ---- Specify model parameters (not required but advised)  -------------------------------------------------
iterations <- 7000      # No of iterations
warmup <- 2000 # Warmup samples to be discarded
init <- 0.5
control <- list(
  adapt_engaged = TRUE,
  adapt_delta = 0.999, #increased from default of 0.8
  stepsize = 0.05, # 0.05 default
  max_treedepth = 15
)

## ---- Set priors -------------------------------------------------
prior <- c(
  # Priors for intercept, and population effects
  prior(normal(0, 0.25), class = "b"), # Sets for factors and their interaction
  prior(normal(0, 1), class = "b", coef = "Intercept"), # Sets for intercept specifically

  # Priors for random effects
  prior(normal(0, 1), class = "sd"),

  # Prior for residual standard deviation
  prior(exponential(2), class = "sigma")
)

## ---- Fit the Bayesian mixed-effects model -------------------------------------------------
# 1. Full Interaction Model (no random intercept)
finalModel_prior <- brm(g_ef.s ~ 0 + Intercept + group * timepoint * condition + (1 + timepoint + group + condition|| subject),
             data = data, 
             family = Gaussian(),  # Specify the appropriate likelihood for your data
             chains = 4,  # Number of chains for MCMC sampling
             cores = 4,   # Number of CPU cores to use
             iter = iterations, # Number of MCMC iterations
             warmup = warmup,
             prior = prior,
             init = init,
             sample_prior = 'only',
             save_pars = save_pars(all = TRUE),
             seed = 22,
             file = "c:/Users/j_m289/Pictures/phd/3. Data Analysis/studies/OKTOS/resting_te/analysis/brms/models/geff/geff_finalModel_prior", # Uncomment to save model object
             control = control # adjust if divergence issues emerge
             )      

finalModel <- brm(g_ef.s ~ 0 + Intercept + group * timepoint * condition + (1 + timepoint + group + condition|| subject),
             data = data, 
             family = Gaussian(),  # Specify the appropriate likelihood for your data
             chains = 4,  # Number of chains for MCMC sampling
             cores = 4,   # Number of CPU cores to use
             iter = iterations, # Number of MCMC iterations
             warmup = warmup,
             prior = prior,
             init = init,
             sample_prior = 'yes',
             save_pars = save_pars(all = TRUE),
             seed = 22,
             file = "c:/Users/j_m289/Pictures/phd/3. Data Analysis/studies/OKTOS/resting_te/analysis/brms/models/geff/geff_finalModel", # Uncomment to save model object
             control = control # adjust if divergence issues emerge
             )                         
finalModel <- readRDS("c:/Users/j_m289/Pictures/phd/3. Data Analysis/studies/OKTOS/resting_te/analysis/brms/models/geff/geff_finalModel.rds")

## ---- Add information criterion to model object  -------------------------------------------------
finalModel <- add_criterion(finalModel, c("waic", "loo"), save_psis =TRUE, moment_match = TRUE, overwrite = TRUE)

## ---- Model Variance -------------------------------------------------
temp1 <- bayes_R2(finalModel) %>%as.data.frame()

loo <- finalModel$criteria$loo$estimates  %>% as.data.frame()
row.names(loo)[row.names(loo) == "model1"] <- finalModel$formula

## ----- Data Export --------------
# Export tables
write.xlsx(temp1, outputFile, sheetName="BayesR2",  append=TRUE)
write.xlsx(loo, outputFile, sheetName="loo",  append=TRUE)