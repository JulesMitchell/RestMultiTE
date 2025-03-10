### MODEL SPECIFICATION AND VALIDATION ###
library(brms)
library(dplyr)
library(xlsx)

outputFile <- "c:/Users/j_m289/Pictures/phd/3. Data Analysis/studies/OKTOS/resting_te/analysis/brms/output/geff_modelcomp.xlsx"

## ----Load Datasets, Clean and Set Factor Levels -------------------------------------------------
#1. Load data
data <- read.csv("c:/Users/j_m289/Pictures/phd/3. Data Analysis/studies/OKTOS/resting_te/analysis/spreadsheets/OKTOS_te_global.csv") %>%
    mutate(
        timepoint = factor(timepoint, levels = c("ses-01", "ses-02", "ses-03"), labels = c("BAS", "POST", "FUP")),
        group = factor(group, levels = c(0, 1), labels = c("NR", "RESP")),
        condition = factor(condition),
        Sex = factor(Sex, levels = c(0, 1), labels = c("Female", "Male"))
      )

# 2. Scaling and plot
# Rescale DV 
model_data$g_ef.s <- scale(model_data$g_ef)[, 1]

### ANALYSIS SECTION ###
## ---- Specify model parameters (not required but advised)  -------------------------------------------------
iterations <- 4000      # No of iterations
warmup <- iterations / 2 # Warmup samples to be discarded
init <- 0.5
control <- list(
  adapt_engaged = TRUE,
  adapt_delta = 0.999, #increased from default of 0.8
  stepsize = 0.05, # 0.05 default # Need to try 0.001
  max_treedepth = 15
)

## ---- Set priors ------------------------------------------------
# For models without random slope correlations
prior <- c(
  # Priors for intercept, and population effects
  prior(normal(0, 0.25), class = "b"), # Sets for factors and their interaction
  prior(normal(0, 1), class = "b", coef = "Intercept"), # Sets for intercept specifically

  # Priors for random effects
  prior(normal(0, 1), class = "sd"),
  prior(lkj(2), class = "cor"),

  # Prior for residual standard deviation
  prior(exponential(2), class = "sigma")
)

prior_cor <- c(
  # Priors for intercept, and population effects
  prior(normal(0, 0.25), class = "b"), # Sets for factors and their interaction
  prior(normal(0, 1), class = "b", coef = "Intercept"), # Sets for intercept specifically

  # Priors for random effects
  prior(normal(0, 1), class = "sd"),

  # Prior for residual standard deviation
  prior(exponential(2), class = "sigma")
)

## ---- Fit the Bayesian mixed-effects model -------------------------------------------------
model1 <- brm(g_ef.s ~ 0 + Intercept + group * timepoint * condition + (1 || subject),
              data = data, 
              family = Gaussian(),  # Specify the appropriate likelihood for your data
              chains = 2,  # Number of chains for MCMC sampling
              cores = 4,   # Number of CPU cores to use
              iter = iterations,  # Number of MCMC iterations
              warmup = warmup,
              prior = prior,
              init = init,
              sample_prior = 'yes',
              save_pars = save_pars(all = TRUE),
              seed = 22,
              file = "c:/Users/j_m289/Pictures/phd/3. Data Analysis/studies/OKTOS/resting_te/analysis/brms/models/geff/geff_model1", # Uncomment save model object
              control = control # adjust if divergence issues emerge
)
 
model2 <- brm(g_ef.s ~ 0 + Intercept + group * timepoint * condition + (1 + timepoint|| subject), 
              data = data, 
              family = Gaussian(),  # Specify the appropriate likelihood for your data
              chains = 2,  # Number of chains for MCMC sampling
              cores = 4,   # Number of CPU cores to use
              iter = iterations,  # Number of MCMC iterations
              warmup = warmup,
              prior = prior,
              init = init,
              sample_prior = 'yes',
              save_pars = save_pars(all = TRUE),
              seed = 22,
              file = "c:/Users/j_m289/Pictures/phd/3. Data Analysis/studies/OKTOS/resting_te/analysis/brms/models/geff/geff_model2", # Uncomment save model object
              control = control # adjust if divergence issues emerge
)

model3 <- brm(g_ef.s ~ 0 + Intercept + group * timepoint * condition + (1 + timepoint + condition|| subject), 
              data = data, 
              family = Gaussian(),  # Specify the appropriate likelihood for your data
              chains = 2,  # Number of chains for MCMC sampling
              cores = 4,   # Number of CPU cores to use
              iter = iterations,  # Number of MCMC iterations
              warmup = warmup,
              prior = prior,
              init = init,
              sample_prior = 'yes',
              save_pars = save_pars(all = TRUE),
              seed = 22,
              file = "c:/Users/j_m289/Pictures/phd/3. Data Analysis/studies/OKTOS/resting_te/analysis/brms/models/geff/geff_model", # Uncomment save model object
              control = control # adjust if divergence issues emerge
)

model4 <- brm(g_ef.s ~ 0 + Intercept + group * timepoint * condition + (1 + timepoint + group || subject), 
              data = data, 
              family = Gaussian(),  # Specify the appropriate likelihood for your data
              chains = 2,  # Number of chains for MCMC sampling
              cores = 4,   # Number of CPU cores to use
              iter = iterations,  # Number of MCMC iterations
              warmup = warmup,
              prior = prior,
              init = init,
              sample_prior = 'yes',
              save_pars = save_pars(all = TRUE),
              seed = 22,
              file = "c:/Users/j_m289/Pictures/phd/3. Data Analysis/studies/OKTOS/resting_te/analysis/brms/models/geff/geff_model4", # Uncomment save model object
              control = control # adjust if divergence issues emerge
)

model5 <- brm(g_ef.s ~ 0 + Intercept + group * timepoint * condition + (1 + timepoint + condition + group|| subject), 
              data = data, 
              family = Gaussian(),  # Specify the appropriate likelihood for your data
              chains = 2,  # Number of chains for MCMC sampling
              cores = 4,   # Number of CPU cores to use
              iter = iterations,  # Number of MCMC iterations
              warmup = warmup,
              prior = prior,
              init = init,
              sample_prior = 'yes',
              save_pars = save_pars(all = TRUE),
              seed = 22,
              file = "c:/Users/j_m289/Pictures/phd/3. Data Analysis/studies/OKTOS/resting_te/analysis/brms/models/geff/geff_model5", # Uncomment save model object
              control = control # adjust if divergence issues emerge
)

## ---- Assess fit (observed vs. predicted) -------------------------------------------------
temp1 <- bayes_R2(model1) %>%as.data.frame()
temp2 <- bayes_R2(model2) %>%as.data.frame()
temp3 <- bayes_R2(model3) %>%as.data.frame()
temp4 <- bayes_R2(model4) %>%as.data.frame()
temp5 <- bayes_R2(model5) %>%as.data.frame()

BayesR2 <- rbind(temp1, temp2, temp3, temp4, temp5)
BayesR2 <- round(BayesR2, 4)

row.names(BayesR2) <- paste("Model", 1:nrow(BayesR2))

## ---- Model comparison  -------------------------------------------------
loo <- loo(model1, model2, model3, model4, model5, compare = TRUE) # requires same numbers of observations
modelcomparison <- loo$diffs %>% as.data.frame()

## ----- Data Export --------------
# Export tables
write.xlsx(BayesR2, outputFile, sheetName="bayesR2",  append=TRUE)
write.xlsx(modelcomparison, outputFile, sheetName="modelselection",  append=TRUE)