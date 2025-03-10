### Model convergence and assumption checks ###
## ----pkgload, include = FALSE-------------------------------------------------
library(ggplot2)
library(brms)
library(tidyverse) # needed for data manipulation
library(dplyr)

metrics <- c('OutDgr') # ClCoef', 'InDgr', 'OutDgr'

for (metric in metrics) {

  # Load model object and extract data
  finalModel <- readRDS(paste0("c:/Users/j_m289/Pictures/phd/3. Data Analysis/studies/OKTOS/resting_te/analysis/brms/models/", tolower(metric), "/", tolower(metric), "_finalmodel.rds.rds"))

  data <- read.csv("c:/Users/j_m289/Pictures/phd/3. Data Analysis/studies/OKTOS/resting_te/analysis/spreadsheets/OKTOS_te_local.csv") %>%
    filter(measure == metric) %>%
    mutate(
      timepoint = factor(timepoint, levels = c("ses-01", "ses-02", "ses-03"), labels = c("BAS", "POST", "FUP")),
      group = factor(group, levels = c(0, 1), labels = c("NR", "RESP")),
      condition = factor(condition),
      channels = relevel(factor(channels), 'Cz')
    )

  data$value.s <- scale(data$value)[, 1]

  # Extract posterior samples for inference
  draws_fit <- as_draws_df(finalModel)

  # Save rhat values
  rhat_vals <- rhat(finalModel)

  # Save Number of effective samples
  neff_vals <- neff_ratio(finalModel)

  ## Plotting
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
    labs(title = paste0(metric, " Rhat Plot")) +
    theme(legend.position.inside = c(.95, .2))

  # 2. Neff
  fig2 <- mcmc_neff_hist(neff_vals) + theme_bw() +
    labs(title = paste0(metric, " Effective Sampling Plot")) +
    theme(legend.position.inside = c(.95, .2))

  # 3. Autocorrelation plot (adjust based on size of model)
  fig3a <- mcmc_acf(draws_fit, 
           pars = vars(b_Intercept:"b_groupRESP:timepointPOST"),
           lags = 5) +
    labs(title = paste0(metric, " Autocorrelation Plots")) +
    theme(legend.position.inside = c(.95, .2))

  fig3b <- mcmc_acf(draws_fit, 
           pars = vars("b_groupRESP:timepointFUP":"b_groupRESP:timepointFUP:conditionEO"),
           lags = 5) +
    labs(title = paste0(metric, " Autocorrelation Plots")) +
    theme(legend.position.inside = c(.95, .2))

  # 4. Plot other diagnostics
  fig4a <- mcmc_trace(draws_fit, 
             pars = vars(b_Intercept:"b_groupRESP:timepointPOST")) +
    labs(title = paste0(metric, " Trace Plots")) +
    theme(legend.position.inside = c(.95, .2))

  fig4b <- mcmc_trace(draws_fit, 
             pars = vars("b_groupRESP:timepointFUP":"b_groupRESP:timepointFUP:conditionEO")) +
    labs(title = paste0(metric, " Trace Plots")) +
    theme(legend.position = c(.95, .2))

  fig5 <- ppc_pit_ecdf(pit=loo_pit(y = data$value.s,
                         yrep = posterior_predict(finalModel),
                         lw = weights(loo(finalModel)$psis_object)))
  fig6 <- pp_check(finalModel, type = "scatter_avg", ndraws = 100)
  fig7 <- pp_check(finalModel, type = "error_hist", ndraws = 12) 

  # 5. Plot posterior distribution checks
  fig8 <- pp_check(finalModel, ndraws = 100)

  # Grouped checks
  fig9 <- pp_check(finalModel, type = "stat_grouped", stat = 'min', group = "timepoint")
  fig10 <- pp_check(finalModel, type = "stat_grouped", stat = 'mean', group = "timepoint")
  fig11 <- pp_check(finalModel, type = "stat_grouped", stat = 'max', group = "timepoint")

  fig12 <- pp_check(finalModel, type = "stat_grouped", stat = 'min', group = "condition")
  fig13 <- pp_check(finalModel, type = "stat_grouped", stat = 'mean', group = "condition")
  fig14 <- pp_check(finalModel, type = "stat_grouped", stat = 'max', group = "condition")

  fig15 <- pp_check(finalModel, type = "stat_grouped", stat = 'min', group = "group")
  fig16 <- pp_check(finalModel, type = "stat_grouped", stat = 'mean', group = "group")
  fig17 <- pp_check(finalModel, type = "stat_grouped", stat = 'max', group = "group")

  ## Figure Export
  # Save figures
  ggsave(paste0("c:/Users/j_m289/Pictures/phd/3. Data Analysis/studies/OKTOS/resting_te/analysis/brms/figures/", tolower(metric), "/", tolower(metric), "_finalmodel_Rhat.jpeg"), plot = fig1, width = 10, height = 10)
  ggsave(paste0("c:/Users/j_m289/Pictures/phd/3. Data Analysis/studies/OKTOS/resting_te/analysis/brms/figures/", tolower(metric), "/", tolower(metric), "_finalmodel_Neff.jpeg"), plot = fig2, width = 10, height = 10)
  ggsave(paste0("c:/Users/j_m289/Pictures/phd/3. Data Analysis/studies/OKTOS/resting_te/analysis/brms/figures/", tolower(metric), "/", tolower(metric), "_finalmodel_AutoCorr_A.jpeg"), plot = fig3a, width = 10, height = 10)
  ggsave(paste0("c:/Users/j_m289/Pictures/phd/3. Data Analysis/studies/OKTOS/resting_te/analysis/brms/figures/", tolower(metric), "/", tolower(metric), "_finalmodel_AutoCorr_B.jpeg"), plot = fig3b, width = 10, height = 10)

  ggsave(paste0("c:/Users/j_m289/Pictures/phd/3. Data Analysis/studies/OKTOS/resting_te/analysis/brms/figures/", tolower(metric), "/", tolower(metric), "_finalmodel_Traces_A.jpeg"), plot = fig4a, width = 10, height = 10)
  ggsave(paste0("c:/Users/j_m289/Pictures/phd/3. Data Analysis/studies/OKTOS/resting_te/analysis/brms/figures/", tolower(metric), "/", tolower(metric), "_finalmodel_Traces_B.jpeg"), plot = fig4b, width = 10, height = 10)

  ggsave(paste0("c:/Users/j_m289/Pictures/phd/3. Data Analysis/studies/OKTOS/resting_te/analysis/brms/figures/", tolower(metric), "/", tolower(metric), "_finalmodel_pitecdf.jpeg"), plot = fig5, width = 10, height = 10)
  ggsave(paste0("c:/Users/j_m289/Pictures/phd/3. Data Analysis/studies/OKTOS/resting_te/analysis/brms/figures/", tolower(metric), "/", tolower(metric), "_finalmodel_ObsVsPred.jpeg"), plot = fig6, width = 10, height = 10)
  ggsave(paste0("c:/Users/j_m289/Pictures/phd/3. Data Analysis/studies/OKTOS/resting_te/analysis/brms/figures/", tolower(metric), "/", tolower(metric), "_finalmodel_ErrorDistribution.jpeg"), plot = fig7, width = 10, height = 10)

  ggsave(paste0("c:/Users/j_m289/Pictures/phd/3. Data Analysis/studies/OKTOS/resting_te/analysis/brms/figures/", tolower(metric), "/", tolower(metric), "_finalmodel_PostCheck_DensOverlay.jpeg"), plot = fig8, width = 10, height = 10)

  ggsave(paste0("c:/Users/j_m289/Pictures/phd/3. Data Analysis/studies/OKTOS/resting_te/analysis/brms/figures/", tolower(metric), "/", tolower(metric), "_finalmodel_PostCheck_Timepoint_min.jpeg"), plot = fig9, width = 10, height = 10)
  ggsave(paste0("c:/Users/j_m289/Pictures/phd/3. Data Analysis/studies/OKTOS/resting_te/analysis/brms/figures/", tolower(metric), "/", tolower(metric), "_finalmodel_PostCheck_Timepoint_mean.jpeg"), plot = fig10, width = 10, height = 10)
  ggsave(paste0("c:/Users/j_m289/Pictures/phd/3. Data Analysis/studies/OKTOS/resting_te/analysis/brms/figures/", tolower(metric), "/", tolower(metric), "_finalmodel_PostCheck_Timepoint_max.jpeg"), plot = fig11, width = 10, height = 10)

  ggsave(paste0("c:/Users/j_m289/Pictures/phd/3. Data Analysis/studies/OKTOS/resting_te/analysis/brms/figures/", tolower(metric), "/", tolower(metric), "_finalmodel_PostCheck_condition_min.jpeg"), plot = fig12, width = 10, height = 10)
  ggsave(paste0("c:/Users/j_m289/Pictures/phd/3. Data Analysis/studies/OKTOS/resting_te/analysis/brms/figures/", tolower(metric), "/", tolower(metric), "_finalmodel_PostCheck_condition_mean.jpeg"), plot = fig13, width = 10, height = 10)
  ggsave(paste0("c:/Users/j_m289/Pictures/phd/3. Data Analysis/studies/OKTOS/resting_te/analysis/brms/figures/", tolower(metric), "/", tolower(metric), "_finalmodel_PostCheck_condition_max.jpeg"), plot = fig14, width = 10, height = 10)

  ggsave(paste0("c:/Users/j_m289/Pictures/phd/3. Data Analysis/studies/OKTOS/resting_te/analysis/brms/figures/", tolower(metric), "/", tolower(metric), "_finalmodel_PostCheck_Responder_min.jpeg"), plot = fig15, width = 10, height = 10)
  ggsave(paste0("c:/Users/j_m289/Pictures/phd/3. Data Analysis/studies/OKTOS/resting_te/analysis/brms/figures/", tolower(metric), "/", tolower(metric), "_finalmodel_PostCheck_Responder_mean.jpeg"), plot = fig16, width = 10, height = 10)
  ggsave(paste0("c:/Users/j_m289/Pictures/phd/3. Data Analysis/studies/OKTOS/resting_te/analysis/brms/figures/", tolower(metric), "/", tolower(metric), "_finalmodel_PostCheck_Responder_max.jpeg"), plot = fig17, width = 10, height = 10)
}