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

metrics <- c('InDgr') #'Btwn', 'InDgr', 'OutDgr'

for (metric in metrics) {

  # Use paste0 to dynamically create filenames based on the metric
  model_summary_output <- paste0("c:/Users/j_m289/Pictures/phd/3. Data Analysis/studies/OKTOS/resting_te/analysis/brms/output/", tolower(metric), "_model_summary.xlsx")
  model_marginal_effects_output <- paste0("c:/Users/j_m289/Pictures/phd/3. Data Analysis/studies/OKTOS/resting_te/analysis/brms/output/", tolower(metric), "_marginaleffects.xlsx")
  response_marginal_diff_output <- paste0("c:/Users/j_m289/Pictures/phd/3. Data Analysis/studies/OKTOS/resting_te/analysis/brms/output/", tolower(metric), "_resp_marginaldiff.csv")
  timepoint_marginal_diff_output <- paste0("c:/Users/j_m289/Pictures/phd/3. Data Analysis/studies/OKTOS/resting_te/analysis/brms/output/", tolower(metric), "_timepoint_marginaldiff.csv")

  ## ---- Load data and model object -----------------
  data <- read.csv("c:/Users/j_m289/Pictures/phd/3. Data Analysis/studies/OKTOS/resting_te/analysis/spreadsheets/OKTOS_te_local.csv") %>%
    filter(measure == metric) %>%
    mutate(
      timepoint = factor(timepoint, levels = c("ses-01", "ses-02", "ses-03"), labels = c("BAS", "POST", "FUP")),
      group = factor(group, levels = c(0, 1), labels = c("NR", "RESP")),
      condition = factor(condition),
      channels = relevel(factor(channels), 'Cz')
    )

  # Load RDS File
  finalModel <- readRDS(paste0("c:/Users/j_m289/Pictures/phd/3. Data Analysis/studies/OKTOS/resting_te/analysis/brms/models/", tolower(metric), "/", tolower(metric), "_finalmodel.rds.rds"))

  # Load function for generating table of comparisons 
  source('c:/Users/j_m289/Pictures/phd/3. Data Analysis/code/functions.R')

  # Get posterior expected values for the data used to create the models. Alternatively, one could use simulated data or a new dataset.
  # Note, epred_draws (i.e., posterior_epred()) used to investigate conditional effects for existing groups (i.e. regions). 
  # This gives the average predicted outcomes for existing groups, incorporating group-specific deviations in intercept or slope. 

  fitted_draws <- finalModel %>%
      epred_draws(newdata = data, re_formula = NULL, ndraws = 1000) %>% # draws limited to 100 for testing
      ungroup() %>%
      select(-.row, -.chain, -.iteration)

  # Extract marginal means
  marginal_means <- fitted_draws %>%
    group_by(group, timepoint, condition, channels) %>%
    summarise(
      mean = mean(.epred),
      HDI_low_0.95 = hdi(.epred, ci=0.95)$CI_low,
      HDI_high_0.95 = hdi(.epred, ci=0.95)$CI_high
    ) 
  marginal_means <- marginal_means %>% as.data.frame()

  ## ---- Generate comparisons from posterior draws -----------
  # 1. Compare channels across responders and non-responders.
  baseline_EC_comp <- fitted_draws %>%
      filter(timepoint == "BAS") %>%
      filter(condition == "EC") %>%
      group_by(group, channels, .draw) %>%
      summarise(
          .epred = mean(.epred),
          .groups = "drop"
      ) %>%
      group_by(channels) %>%
      compare_levels(variable = .epred, by = group, comparison = "pairwise") %>%
      calc_compare_levels_stats(., .epred)

  baseline_EO_comp <- fitted_draws %>%
      filter(timepoint == "BAS") %>%
      filter(condition == "EO") %>%
      group_by(group, channels, .draw) %>%
      summarise(
          .epred = mean(.epred),
          .groups = "drop"
      ) %>%
      group_by(channels) %>%
      compare_levels(variable = .epred, by = group, comparison = "pairwise") %>%
      calc_compare_levels_stats(., .epred)

  # 1b. Post-treatment
  post_EC_comp <- fitted_draws %>%
      filter(timepoint == "POST") %>%
      filter(condition == "EC") %>%
      group_by(group, channels, .draw) %>%
      summarise(
          .epred = mean(.epred),
          .groups = "drop"
      ) %>%
      group_by(channels) %>%
      compare_levels(variable = .epred, by = group, comparison = "pairwise") %>%
      calc_compare_levels_stats(., .epred)

  post_EO_comp <- fitted_draws %>%
      filter(timepoint == "POST") %>%
      filter(condition == "EO") %>%
      group_by(group, channels, .draw) %>%
      summarise(
          .epred = mean(.epred),
          .groups = "drop"
      ) %>%
      group_by(channels) %>%
      compare_levels(variable = .epred, by = group, comparison = "pairwise") %>%
      calc_compare_levels_stats(., .epred)	

  # 1c. Follow-up
  fup_EC_comp <- fitted_draws %>%
      filter(timepoint == "FUP") %>%
      filter(condition == "EC") %>%
      group_by(group, channels, .draw) %>%
      summarise(
          .epred = mean(.epred),
          .groups = "drop"
      ) %>%
      group_by(channels) %>%
      compare_levels(variable = .epred, by = group, comparison = "pairwise") %>%
      calc_compare_levels_stats(., .epred)

  fup_EO_comp <- fitted_draws %>%
      filter(timepoint == "FUP") %>%
      filter(condition == "EO") %>%
      group_by(group, channels, .draw) %>%
      summarise(
          .epred = mean(.epred),
          .groups = "drop"
      ) %>%
      group_by(channels) %>%
      compare_levels(variable = .epred, by = group, comparison = "pairwise") %>%
      calc_compare_levels_stats(., .epred)

  # 2. Compare timepoint for eyes closed across responders and non-responders.
  timepoint_RESP_EC_comp <- fitted_draws %>%
    filter(group != "NR") %>%
    filter(condition == "EC") %>%
      group_by(timepoint, channels, .draw) %>%
      summarise(
          .epred = mean(.epred),
          .groups = "drop"
      ) %>%
      group_by(channels) %>%
      compare_levels(variable = .epred, by = timepoint, comparison = "pairwise") %>%
      calc_compare_levels_stats(., .epred)

  timepoint_RESP_EO_comp <- fitted_draws %>%
    filter(group != "NR") %>%
    filter(condition == "EO") %>%
      group_by(timepoint, channels, .draw) %>%
      summarise(
          .epred = mean(.epred),
          .groups = "drop"
      ) %>%
      group_by(channels) %>%
      compare_levels(variable = .epred, by = timepoint, comparison = "pairwise") %>%
      calc_compare_levels_stats(., .epred)    

  timepoint_NR_EC_comp <- fitted_draws %>%
    filter(group != "RESP") %>%
    filter(condition == "EC") %>%
      group_by(timepoint, channels, .draw) %>%
      summarise(
          .epred = mean(.epred),
          .groups = "drop"
      ) %>%
      group_by(channels) %>%
      compare_levels(variable = .epred, by = timepoint, comparison = "pairwise") %>%
      calc_compare_levels_stats(., .epred)

  timepoint_NR_EO_comp <- fitted_draws %>%
    filter(group != "RESP") %>%
    filter(condition == "EO") %>%
      group_by(timepoint, channels, .draw) %>%
      summarise(
          .epred = mean(.epred),
          .groups = "drop"
      ) %>%
      group_by(channels) %>%
      compare_levels(variable = .epred, by = timepoint, comparison = "pairwise") %>%
      calc_compare_levels_stats(., .epred) 

  ## ----- Data Export --------------
  # Export tables
  # write.xlsx(effects, model_summary_output, sheetName="EffectEstimates",  append=TRUE)
  # write.xlsx(pvalues, model_summary_output, sheetName="Pvalues",  append=TRUE)
  write.xlsx(marginal_means, model_marginal_effects_output, sheetName="MarginalMeans",  append=TRUE)

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
      group = 'NR',
      comparison = paste(condition, timepoint, sep = ":"), 
      pvalue = (1 - pd)*2
    )

  timepoint_NR_EO_comp <- timepoint_NR_EO_comp %>%
    mutate(
      condition = "EO",
      group = 'NR',
      comparison = paste(condition, timepoint, sep = ":"), 
      pvalue = (1 - pd)*2
    )

  timepoint_RESP_EC_comp <- timepoint_RESP_EC_comp %>%
    mutate(
      condition = "EC",
      group = 'RESP', 
      comparison = paste(condition, timepoint, sep = ":"), 
      pvalue = (1 - pd)*2
    )

  timepoint_RESP_EO_comp <- timepoint_RESP_EO_comp %>%
    mutate(
      condition = "EO",
      group = 'RESP', 
      comparison = paste(condition, timepoint, sep = ":"), 
      pvalue = (1 - pd)*2
    )

  combined_timepoint_comp <- bind_rows(timepoint_NR_EC_comp, timepoint_NR_EO_comp, timepoint_RESP_EC_comp, timepoint_RESP_EO_comp)
  write.csv(combined_timepoint_comp, timepoint_marginal_diff_output)

  # ## ----- Plotting (adjust code to plot desired combinations) ----------
  # fitted_draws_filtered <- fitted_draws %>%
  #   filter(condition == "EC" & timepoint == "BAS") # Subset fitted_draws as desired

  # # Create a new variable to categorize channels (Left, Central, Right)
  # fitted_draws_filtered <- fitted_draws_filtered %>%
  #   mutate(channel_group = case_when(
  #     grepl("\\d$", channels) & as.numeric(substr(channels, nchar(channels), nchar(channels))) %% 2 == 1 ~ "Left",  # Odd numbers
  #     grepl("z$", channels) ~ "Central",  # Ends with 'z'
  #     grepl("\\d$", channels) & as.numeric(substr(channels, nchar(channels), nchar(channels))) %% 2 == 0 ~ "Right",  # Even numbers
  #     TRUE ~ "Other"  # Fallback case
  #   ))

  # # Reorder the channel_group factor so that "Central" is in the middle
  # fitted_draws_filtered$channel_group <- factor(fitted_draws_filtered$channel_group, 
  #                                               levels = c("Left", "Central", "Right"))

  # # Plot using facet_wrap() to separate by channel group (Left, Central, Right)
  # fig1 <- ggplot(fitted_draws_filtered, aes(x = .epred, y = channels, fill = group)) +
  #         stat_halfeye() +
  #         labs(x = "Betweeness Centrality", y = NULL,
  #         fill = "Response Status",
  #         subtitle = "Posterior predictions by Channel Group") +
  #         facet_wrap(~ channel_group, scales = "free_y", ncol = 3) +  # Split into 3 columns (Left, Central, Right)
  #         theme(legend.position = "bottom", 
  #         strip.text = element_text(size = 8),  # Smaller facet labels
  #         axis.text = element_text(size = 8),   # Smaller axis labels
  #         panel.spacing = unit(0.5, "lines"))   # Reduce spacing between panels

  # # Save figures
  # ggsave("c:/Users/j_m289/Pictures/phd/3. Data Analysis/studies/OKTOS/resting_te/analysis/brms/figures/", metric, "_Interaction1.jpeg", plot = fig1, width = 10, height = 10)
}

