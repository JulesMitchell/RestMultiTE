### PRIOR PREDICTIVE CHECK ###
## ----pkgload, include = FALSE-------------------------------------------------
library(brms)
library(ggplot2)

metrics <- c('Btwn', 'ClCoef', 'InDgr', 'OutDgr')

for (metric in metrics) {

	## ---- Load saved model object -----------------
	finalModel_prior <- readRDS(paste0("c:/Users/j_m289/Pictures/phd/3. Data Analysis/studies/OKTOS/resting_te/analysis/brms/models/", tolower(metric), "/", metric, "_finalmodel_prior.rds"))

	## ---- Review model prior and summary(estimates and CIs) -------------------------------------------------
	#. 1 Prior predictive checks
	prior_summary(finalModel_prior)  # Re-specify priors and re-run model if needed
	summary(finalModel_prior) # Check estimates produced from priors

	# if MV use resp argument, else remove. 
	fig1 <- pp_check(finalModel_prior, type  = 'hist', ndraws = 12, prefix='ppd')

	## ---- Export figures -------------------------------------------------
	ggsave(paste0("c:/Users/j_m289/Pictures/phd/3. Data Analysis/studies/OKTOS/resting_te/analysis/brms/figures/", tolower(metric), "/", tolower(metric), "_prior_PPCheck_hist.jpeg"), plot = fig1, width = 10, height = 10)
}