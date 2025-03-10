### PRIOR PREDICTIVE CHECK ###
## ----pkgload, include = FALSE-------------------------------------------------
library(brms)
library(ggplot2)

## ---- Load saved model object (useful for coming back to analysis later) -----------------
# 1. To load
finalModel_prior <- readRDS("c:/Users/j_m289/Pictures/phd/3. Data Analysis/studies/OKTOS/resting_te/analysis/brms/models/geff/geff_finalModel_prior.rds") # adjust filename as required

## ---- Review model prior and summary(estimates and CIs) -------------------------------------------------
prior_summary(finalModel_prior)  # Re-specify priors and re-run model if needed
summary(finalModel_prior) # Check estimates produced from priors

fig1 <- pp_check(finalModel_prior, type  = 'hist', ndraws = 12, prefix='ppd')

## ---- Export figures -------------------------------------------------
ggsave("c:/Users/j_m289/Pictures/phd/3. Data Analysis/studies/OKTOS/resting_te/analysis/brms/figures/geff/geff_prior_PPCheck_hist.jpeg", plot = fig1, width = 10, height = 10)