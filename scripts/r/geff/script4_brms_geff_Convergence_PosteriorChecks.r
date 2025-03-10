### Model convergence and assumption checks ###
## ----pkgload, include = FALSE-------------------------------------------------
library(ggplot2)
library(brms)
library(tidyverse) # needed for data manipulation

## ---- Load model object and save rhat and neff values -----------------
# 1. To load
finalModel <- readRDS("c:/Users/j_m289/Pictures/phd/3. Data Analysis/studies/OKTOS/resting_te/analysis/brms/models/geff/geff_finalModel.rds") # adjust filename as required

data <- read.csv("c:/Users/j_m289/Pictures/phd/3. Data Analysis/studies/OKTOS/resting_te/analysis/spreadsheets/OKTOS_te_global.csv") %>%
    mutate(
        timepoint = factor(timepoint, levels = c("ses-01", "ses-02", "ses-03"), labels = c("BAS", "POST", "FUP")),
        group = factor(group, levels = c(0, 1), labels = c("NR", "RESP")),
        condition = factor(condition))

# 2. Scaling and plot
# Rescale DV 
data$g_ef.s <- scale(data$g_ef)[, 1]

# 2. Extract posterior samples for inference
draws_fit <- as_draws_df(finalModel)

# 3. Save rhat values
rhat_vals <- rhat(finalModel)

# 4. Save Number of effective samples
neff_vals <- neff_ratio(finalModel)

## ---- Plotting-----------------
library(bayesplot) # Loaded after rhat as it is masked by BRMS

# Define function for pit_ecdf plotting
loo_pit <- function(y, yrep, lw) {
  pit <- vapply(seq_len(ncol(yrep)), function(j) {
    sel_min <- yrep[, j] < y[j]
    pit_min <- exp_log_sum_exp(lw[sel_min,j])
    sel_sup <- yrep[, j] == y[j]
    pit_sup <- pit_min + exp_log_sum_exp(lw[sel_sup,j])
    runif(1, pit_min, pit_sup)
  }, FUN.VALUE = 1)
  pmax(pmin(pit, 1), 0)
}
exp_log_sum_exp <- function(x) {
  m <- suppressWarnings(max(x))
  exp(m + log(sum(exp(x - m))))
}

# 1. rhat  
fig1 <- mcmc_rhat_hist(rhat_vals) + theme_bw() +
  labs(title = "Rhat Plot") +
  theme(legend.position.inside = c(.95, .2))

# 2. Neff
fig2 <- mcmc_neff_hist(neff_vals)  + theme_bw() +
  labs(title = "Effective Sampling Plot") +
  theme(legend.position.inside = c(.95, .2))

# 3. Autocorrelation plot (adjust based on size of model)
fig3a <- mcmc_acf(draws_fit, 
         pars = vars(b_Intercept:"b_groupRESP:timepointPOST"),
         lags = 5) +
  labs(title = "Autocorrelation Plots") +
  theme(legend.position.inside = c(.95, .2))

fig3b <- mcmc_acf(draws_fit, 
         pars = vars("b_groupRESP:timepointFUP":"b_groupRESP:timepointFUP:conditionEO"),
         lags = 5) +
  labs(title = "Autocorrelation Plots") +
  theme(legend.position.inside = c(.95, .2))

# 4. Plot other diagnostics
fig4a <- mcmc_trace(draws_fit, 
           pars = vars(b_Intercept:"b_groupRESP:timepointPOST")) +
  labs(title = "Trace Plots") +
  theme(legend.position.inside = c(.95, .2))

fig4b <- mcmc_trace(draws_fit, 
           pars = vars("b_groupRESP:timepointFUP":"b_groupRESP:timepointFUP:conditionEO")) +
  labs(title = "Trace Plots") +
  theme(legend.position.inside = c(.95, .2))

fig5 <- ppc_pit_ecdf(pit=loo_pit(y = data$g_ef.s,
                         yrep = posterior_predict(finalModel),
                         lw = weights(loo(finalModel)$psis_object)))
fig6 <- pp_check(finalModel, type = "scatter_avg", ndraws = 100)
fig7 <- pp_check(finalModel, type = "error_hist", ndraws = 12) 

# 5. Plot posterior distribution checks
# Overall
fig8 <- pp_check(finalModel, ndraws = 100)

# Grouped checks
fig9 <- pp_check(finalModel, type = "stat_grouped", stat = 'min', group = "timepoint") # Keep (use for priors and posteriors)
fig10 <- pp_check(finalModel, type = "stat_grouped", stat = 'mean', group = "timepoint") # Keep (use for priors and posteriors)
fig11 <- pp_check(finalModel, type = "stat_grouped", stat = 'max', group = "timepoint") # Keep (use for priors and posteriors)

fig12 <- pp_check(finalModel, type = "stat_grouped", stat = 'min', group = "group") # Keep (use for priors and posteriors)
fig13 <- pp_check(finalModel, type = "stat_grouped", stat = 'mean', group = "group") # Keep (use for priors and posteriors)
fig14 <- pp_check(finalModel, type = "stat_grouped", stat = 'max', group = "group") # Keep (use for priors and posteriors)

fig15 <- pp_check(finalModel, type = "stat_grouped", stat = 'min', group = "condition") # Keep (use for priors and posteriors)
fig16 <- pp_check(finalModel, type = "stat_grouped", stat = 'mean', group = "condition") # Keep (use for priors and posteriors)
fig17 <- pp_check(finalModel, type = "stat_grouped", stat = 'max', group = "condition") # Keep (use for priors and posteriors)

## ---- Figure Export -----------------
# Save figures
ggsave("c:/Users/j_m289/Pictures/phd/3. Data Analysis/studies/OKTOS/resting_te/analysis/brms/figures/geff/geff_finalmodel_Rhat.jpeg", plot = fig1, width = 10, height = 10)
ggsave("c:/Users/j_m289/Pictures/phd/3. Data Analysis/studies/OKTOS/resting_te/analysis/brms/figures/geff/geff_finalmodel_Neff.jpeg", plot = fig2, width = 10, height = 10)
ggsave("c:/Users/j_m289/Pictures/phd/3. Data Analysis/studies/OKTOS/resting_te/analysis/brms/figures/geff/geff_finalmodel_AutoCorr_A.jpeg", plot = fig3a, width = 10, height = 10)
ggsave("c:/Users/j_m289/Pictures/phd/3. Data Analysis/studies/OKTOS/resting_te/analysis/brms/figures/geff/geff_finalmodel_AutoCorr_B.jpeg", plot = fig3b, width = 10, height = 10)

ggsave("c:/Users/j_m289/Pictures/phd/3. Data Analysis/studies/OKTOS/resting_te/analysis/brms/figures/geff/geff_finalmodel_Traces_A.jpeg", plot = fig4a, width = 10, height = 10)
ggsave("c:/Users/j_m289/Pictures/phd/3. Data Analysis/studies/OKTOS/resting_te/analysis/brms/figures/geff/geff_finalmodel_Traces_B.jpeg", plot = fig4b, width = 10, height = 10)

ggsave("c:/Users/j_m289/Pictures/phd/3. Data Analysis/studies/OKTOS/resting_te/analysis/brms/figures/geff/geff_finalmodel_pit_ecdf.jpeg", plot = fig5, width = 10, height = 10)
ggsave("c:/Users/j_m289/Pictures/phd/3. Data Analysis/studies/OKTOS/resting_te/analysis/brms/figures/geff/geff_finalmodel_ObsVsPred.jpeg", plot = fig6, width = 10, height = 10)
ggsave("c:/Users/j_m289/Pictures/phd/3. Data Analysis/studies/OKTOS/resting_te/analysis/brms/figures/geff/geff_finalmodel_ErrorDistribution.jpeg", plot = fig7, width = 10, height = 10)

ggsave("c:/Users/j_m289/Pictures/phd/3. Data Analysis/studies/OKTOS/resting_te/analysis/brms/figures/geff/geff_finalmodel_PostCheck_DensOverlay.jpeg", plot = fig8, width = 10, height = 10)

ggsave("c:/Users/j_m289/Pictures/phd/3. Data Analysis/studies/OKTOS/resting_te/analysis/brms/figures/geff/geff_finalmodel_PostCheck_Timepoint_min.jpeg", plot = fig9, width = 10, height = 10)
ggsave("c:/Users/j_m289/Pictures/phd/3. Data Analysis/studies/OKTOS/resting_te/analysis/brms/figures/geff/geff_finalmodel_PostCheck_Timepoint_mean.jpeg", plot = fig10, width = 10, height = 10)
ggsave("c:/Users/j_m289/Pictures/phd/3. Data Analysis/studies/OKTOS/resting_te/analysis/brms/figures/geff/geff_finalmodel_PostCheck_Timepoint_max.jpeg", plot = fig11, width = 10, height = 10)

ggsave("c:/Users/j_m289/Pictures/phd/3. Data Analysis/studies/OKTOS/resting_te/analysis/brms/figures/geff/geff_finalmodel_PostCheck_Responder_min.jpeg", plot = fig12, width = 10, height = 10)
ggsave("c:/Users/j_m289/Pictures/phd/3. Data Analysis/studies/OKTOS/resting_te/analysis/brms/figures/geff/geff_finalmodel_PostCheck_Responder_mean.jpeg", plot = fig13, width = 10, height = 10)
ggsave("c:/Users/j_m289/Pictures/phd/3. Data Analysis/studies/OKTOS/resting_te/analysis/brms/figures/geff/geff_finalmodel_PostCheck_Responder_max.jpeg", plot = fig14, width = 10, height = 10)

ggsave("c:/Users/j_m289/Pictures/phd/3. Data Analysis/studies/OKTOS/resting_te/analysis/brms/figures/geff/geff_finalmodel_PostCheck_Condition_min.jpeg", plot = fig15, width = 10, height = 10)
ggsave("c:/Users/j_m289/Pictures/phd/3. Data Analysis/studies/OKTOS/resting_te/analysis/brms/figures/geff/geff_finalmodel_PostCheck_Condition_mean.jpeg", plot = fig16, width = 10, height = 10)
ggsave("c:/Users/j_m289/Pictures/phd/3. Data Analysis/studies/OKTOS/resting_te/analysis/brms/figures/geff/geff_finalmodel_PostCheck_Condition_max.jpeg", plot = fig17, width = 10, height = 10)