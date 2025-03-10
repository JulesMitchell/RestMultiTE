### SUMMARRY AND POST HOC CONTRASTS###
## ----pkgload, include = FALSE-------------------------------------------------
library(ggplot2)
library(brms)
library(tidyr)
library(tidyverse) # needed for data manipulation
library(tidybayes) # needed for data manipulation
library(emmeans)
library(bayestestR)
library(knitr)
library(xlsx)

## ---- Set file location and name for data export -------------------------------------------------
model_summary_output <- "c:/Users/j_m289/Pictures/phd/3. Data Analysis/studies/OKTOS/resting_te/analysis/brms/output/geff_model_summary.xlsx"
model_marginal_effects_output <- "c:/Users/j_m289/Pictures/phd/3. Data Analysis/studies/OKTOS/resting_te/analysis/brms/output/geff_marginaleffects.xlsx"
response_marginal_diff_output <- "c:/Users/j_m289/Pictures/phd/3. Data Analysis/studies/OKTOS/resting_te/analysis/brms/output/geff_resp_marginaldiff.csv"
timepoint_marginal_diff_output <- "c:/Users/j_m289/Pictures/phd/3. Data Analysis/studies/OKTOS/resting_te/analysis/brms/output/geff_timepoint_marginaldiff.csv"

## ---- Load dataframe and model object. Define functions and plotting theme. 
# Load data and filter 
data <- read.csv("c:/Users/j_m289/Pictures/phd/3. Data Analysis/studies/OKTOS/resting_te/analysis/spreadsheets/OKTOS_te_global.csv") %>%
  mutate(
    timepoint = factor(timepoint, levels = c("ses-01", "ses-02", "ses-03"), labels = c("BAS", "POST", "FUP")),
    group = factor(group, levels = c(0, 1), labels = c("NR", "RESP")),
    condition = factor(condition)
  )

# Load RDS File
finalModel <- readRDS("c:/Users/j_m289/Pictures/phd/3. Data Analysis/studies/OKTOS/resting_te/analysis/brms/models/geff/geff_finalModel.rds")

# Load function for generating table of comparisons 
source('c:/Users/j_m289/Pictures/phd/3. Data Analysis/code/functions.R')

## ----- Extract posterior draws -----------------
# Get posterior expected values for the data used to create the models. Alternatively, one could use simulated data or a new dataset.
# Note, epred_draws (i.e., posterior_epred()) used to investigate conditional effects for existing groups (i.e. regions). 
# This gives the average predicted outcomes for existing groups, incorporating group-specific deviations in intercept or slope. 

fitted_draws <- finalModel %>%
    epred_draws(newdata = data, re_formula = NULL, ndraws = 1000) %>% # draws limited to 100 for testing
    ungroup() %>%
    select(-.row, -.chain, -.iteration)

# Extract marginal means
marginal_means <- fitted_draws %>%
  group_by(group, timepoint, condition) %>%
  summarise(
    mean = mean(.epred),
    HDI_low_0.95 = hdi(.epred, ci=0.95)$CI_low,
    HDI_high_0.95 = hdi(.epred, ci=0.95)$CI_high
  ) 
marginal_means <- marginal_means %>% as.data.frame()

# Plot Marginal Means and Confidence Intervals
pd <- position_dodge(0.5)  # Set the dodge width to move points and error bars horizontally

fig1 <- marginal_means %>%
  ggplot(aes(x = condition, y = mean, color = group)) +
  # Points for the estimated values
  geom_point(size = 3, position = pd) +
  # Error bars for the confidence intervals
  geom_errorbar(aes(ymin = HDI_low_0.95, ymax = HDI_high_0.95), width = 0.2, position = pd) +
  # Labels and title
  labs(x = "Condition", y = "Global Efficiency", color = "Group") +  # Change color legend title to 'Group'
  # Remove grid lines
  theme_light() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 14),  # Increase x-axis text size
    axis.text.y = element_text(size = 14),  # Increase y-axis text size
    axis.title.x = element_text(size = 16),  # Increase x-axis label size
    axis.title.y = element_text(size = 16),  # Increase y-axis label size
    strip.text = element_text(color = "black", size = 12),
    panel.grid = element_blank()  # Removes all grid lines
  ) +
  facet_grid(~ timepoint)


fig2 <- marginal_means %>%
  ggplot(aes(x = timepoint, y = mean, color = group, group = group)) +
  # Points for the estimated values
  geom_point(size = 3, position = pd) +
  # Error bars for the confidence intervals
  geom_errorbar(aes(ymin = HDI_low_0.95, ymax = HDI_high_0.95), width = 0.2, position = pd) +
  # Lines connecting the points
  geom_line(linewidth = 1, position = pd) +
  # Labels and title
  labs(x = "Timepoint", y = "Global Efficiency", color = "Group") +  # Change color legend title to 'Group'
  # Remove grid lines
  theme_light() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 14),  # Increase x-axis text size
    axis.text.y = element_text(size = 14),  # Increase y-axis text size
    axis.title.x = element_text(size = 16),  # Increase x-axis label size
    axis.title.y = element_text(size = 16),  # Increase y-axis label size
    strip.text = element_text(color = "black", size = 12),
    panel.grid = element_blank()  # Removes all grid lines
  ) +
  facet_grid(~ condition)

## ---- Generate median differences from posterior draws -----------
# 1. Compare channels across responders and non-responders.
# 1a. Baseline
baseline_EC_comp <- fitted_draws %>%
    filter(timepoint == "BAS") %>%
	  filter(condition == "EC") %>%
    group_by(group, .draw) %>%
    summarise(
        .epred = mean(.epred),
        .groups = "drop"
    ) %>%
    compare_levels(variable = .epred, by = group, comparison = "pairwise") %>%
    calc_compare_levels_stats(., .epred)

baseline_EO_comp <- fitted_draws %>%
    filter(timepoint == "BAS") %>%
	  filter(condition == "EO") %>%
    group_by(group, .draw) %>%
    summarise(
        .epred = mean(.epred),
        .groups = "drop"
    ) %>%
    compare_levels(variable = .epred, by = group, comparison = "pairwise") %>%
    calc_compare_levels_stats(., .epred)

# 1b. Post-treatment
post_EC_comp <- fitted_draws %>%
    filter(timepoint == "POST") %>%
	  filter(condition == "EC") %>%
    group_by(group, .draw) %>%
    summarise(
        .epred = mean(.epred),
        .groups = "drop"
    ) %>%
    compare_levels(variable = .epred, by = group, comparison = "pairwise") %>%
    calc_compare_levels_stats(., .epred)

post_EO_comp <- fitted_draws %>%
    filter(timepoint == "POST") %>%
	  filter(condition == "EO") %>%
    group_by(group, .draw) %>%
    summarise(
        .epred = mean(.epred),
        .groups = "drop"
    ) %>%
    compare_levels(variable = .epred, by = group, comparison = "pairwise") %>%
    calc_compare_levels_stats(., .epred)	

# 1c. Follow-up
fup_EC_comp <- fitted_draws %>%
    filter(timepoint == "FUP") %>%
	  filter(condition == "EC") %>%
    group_by(group, .draw) %>%
    summarise(
        .epred = mean(.epred),
        .groups = "drop"
    ) %>%
    compare_levels(variable = .epred, by = group, comparison = "pairwise") %>%
    calc_compare_levels_stats(., .epred)

fup_EO_comp <- fitted_draws %>%
    filter(timepoint == "FUP") %>%
	  filter(condition == "EO") %>%
    group_by(group, .draw) %>%
    summarise(
        .epred = mean(.epred),
        .groups = "drop"
    ) %>%
    compare_levels(variable = .epred, by = group, comparison = "pairwise") %>%
    calc_compare_levels_stats(., .epred)

# 2. Compare timepoint for eyes closed across responders and non-responders.
timepoint_RESP_EC_comp <- fitted_draws %>%
	filter(group != "NR") %>%
	filter(condition == "EC") %>%
    group_by(timepoint, .draw) %>%
    summarise(
        .epred = mean(.epred),
        .groups = "drop"
    ) %>%
    compare_levels(variable = .epred, by = timepoint, comparison = "pairwise") %>%
    calc_compare_levels_stats(., .epred)

timepoint_RESP_EO_comp <- fitted_draws %>%
	filter(group != "NR") %>%
	filter(condition == "EO") %>%
    group_by(timepoint, .draw) %>%
    summarise(
        .epred = mean(.epred),
        .groups = "drop"
    ) %>%
    compare_levels(variable = .epred, by = timepoint, comparison = "pairwise") %>%
    calc_compare_levels_stats(., .epred)    

timepoint_NR_EC_comp <- fitted_draws %>%
	filter(group != "RESP") %>%
	filter(condition == "EC") %>%
    group_by(timepoint, .draw) %>%
    summarise(
        .epred = mean(.epred),
        .groups = "drop"
    ) %>%
    compare_levels(variable = .epred, by = timepoint, comparison = "pairwise") %>%
    calc_compare_levels_stats(., .epred)

timepoint_NR_EO_comp <- fitted_draws %>%
	filter(group != "RESP") %>%
	filter(condition == "EO") %>%
    group_by(timepoint, .draw) %>%
    summarise(
        .epred = mean(.epred),
        .groups = "drop"
    ) %>%
    compare_levels(variable = .epred, by = timepoint, comparison = "pairwise") %>%
    calc_compare_levels_stats(., .epred) 

## ----- Data Export --------------
# Export tables
write.xlsx(effects, model_summary_output, sheetName="EffectEstimates",  append=TRUE)
write.xlsx(pvalues, model_summary_output, sheetName="Pvalues",  append=TRUE)
write.xlsx(marginal_means, model_marginal_effects_output)

# Export posterior predictions
# Add a 'condition' column to each table to identify whether it's "EC" or "EO" and the timepoint comparison. 
baseline_EC_comp <- baseline_EC_comp %>%
  mutate(
    condition = "EC", 
	comparison = "BAS: RESP - NR", 
    pvalue = (1 - pd)*2
  )

baseline_EO_comp <- baseline_EO_comp %>%
  mutate(
    condition = "EO", 
	comparison = "BAS: RESP - NR", 
    pvalue = (1 - pd)*2
  )

post_EC_comp <- post_EC_comp %>%
  mutate(
    condition = "EC", 
	comparison = "POST: RESP - NR", 
    pvalue = (1 - pd)*2
  )

post_EO_comp <- post_EO_comp %>%
  mutate(
    condition = "EO", 
	comparison = "POST: RESP - NR", 
    pvalue = (1 - pd)*2
  )

fup_EC_comp <- fup_EC_comp %>%
  mutate(
    condition = "EC", 
	comparison = "FUP: RESP - NR", 
    pvalue = (1 - pd)*2
  )

fup_EO_comp <- fup_EO_comp %>%
  mutate(
    condition = "EO", 
	comparison = "FUP: RESP - NR", 
    pvalue = (1 - pd)*2
  )

# Combine the two tables using bind_rows and save
combined_response_comp <- bind_rows(baseline_EC_comp, baseline_EO_comp, post_EC_comp, post_EO_comp, fup_EC_comp, fup_EO_comp)
write.csv(combined_response_comp, response_marginal_diff_output)

# Now for the timepoint comparisons. 
timepoint_NR_EC_comp <- timepoint_NR_EC_comp %>%
  mutate(
    condition = "EC", 
	comparison = "NR: EC", 
    pvalue = (1 - pd)*2
  )

timepoint_NR_EO_comp <- timepoint_NR_EO_comp %>%
  mutate(
    condition = "EO", 
	comparison = "NR:EO", 
    pvalue = (1 - pd)*2
  )

timepoint_RESP_EC_comp <- timepoint_RESP_EC_comp %>%
  mutate(
    condition = "EC", 
	comparison = "RESP: EC", 
    pvalue = (1 - pd)*2
  )

timepoint_RESP_EO_comp <- timepoint_RESP_EO_comp %>%
  mutate(
    condition = "EO", 
	comparison = "RESP:EO", 
    pvalue = (1 - pd)*2
  )

combined_timepoint_comp <- bind_rows(timepoint_NR_EC_comp, timepoint_NR_EO_comp, timepoint_RESP_EC_comp, timepoint_RESP_EO_comp)
write.csv(combined_timepoint_comp, timepoint_marginal_diff_output)

# Save figures
ggsave("c:/Users/j_m289/Pictures/phd/3. Data Analysis/studies/OKTOS/resting_te/analysis/brms/figures/geff/geff_Interaction1.jpeg", plot = fig1, width = 10, height = 10)
ggsave("c:/Users/j_m289/Pictures/phd/3. Data Analysis/studies/OKTOS/resting_te/analysis/brms/figures/geff/geff_Interaction2.jpeg", plot = fig2, width = 10, height = 10)