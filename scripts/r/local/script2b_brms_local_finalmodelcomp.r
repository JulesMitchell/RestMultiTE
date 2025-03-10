### MODEL Comparison ###
## ----pkgload, include = FALSE-------------------------------------------------
library(brms)
library(dplyr)
library(xlsx)

metrics <- c('Btwn', 'Clcoef', 'OutDgr', 'InDgr') 

for (metric in metrics) {
  
  # Checkpoint 1: Start processing the metric
  print(paste("Processing metric:", metric))
  
  outputFile <- paste0("c:/Users/j_m289/Pictures/phd/3. Data Analysis/studies/OKTOS/resting_te/analysis/brms/output/", tolower(metric), "_finaloutput.xlsx")
  
  ## ----Load RDS files -------------------------------------------------
  print(paste("Loading final model for metric:", metric))
  finalModel <- readRDS(paste0("c:/Users/j_m289/Pictures/phd/3. Data Analysis/studies/OKTOS/resting_te/analysis/brms/models/", tolower(metric), "/", metric, "_finalmodel.rds.rds"))
  
  print(finalModel)

  print(loo(finalModel))
  # # Checkpoint 2: Model loaded
  # print(paste("Model loaded for metric:", metric))
  
  # ## ---- Add information criterion to model object -------------------------------------------------
  # print(paste("Adding information criteria (WAIC and LOO) for metric:", metric))
  # finalModel <- add_criterion(finalModel, c("waic", "loo"), save_psis = TRUE, overwrite = TRUE, file = paste0("c:/Users/j_m289/Pictures/phd/3. Data Analysis/studies/OKTOS/resting_te/analysis/brms/models/", tolower(metric), "/", metric, "_finalmodel.rds"))
  
  # # Checkpoint 3: Information criteria added
  # print(paste("Information criteria added for metric:", metric))
  
  # ## ---- Model Variance -------------------------------------------------
  # print(paste("Calculating Bayes R2 for metric:", metric))
  # temp1 <- bayes_R2(finalModel) %>% as.data.frame()
  
  # loo <- finalModel$criteria$loo$estimates %>% as.data.frame()
  # row.names(loo)[row.names(loo) == "model1"] <- finalModel$formula
  
  # print(paste("Bayes R2 and LOO calculated for metric:", metric))
  
  # print(loo(finalModel))
  
  # ## ----- Data Export --------------
  # print(paste("Exporting data for metric:", metric))
  # # Export tables
  # write.xlsx(temp1, outputFile, sheetName=paste0(metric, "bayesR2"), append=TRUE)
  # write.xlsx(loo, outputFile, sheetName=paste0(metric, "loo"), append=TRUE)
  
  # # Checkpoint 4: Data exported
  # print(paste("Data exported for metric:", metric))

  gc()
  
}

# Final checkpoint after the loop finishes
print("All metrics processed and data exported.")