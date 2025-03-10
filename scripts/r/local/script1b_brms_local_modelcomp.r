### MODEL Comparison ###
## ----pkgload, include = FALSE-------------------------------------------------
library(brms)
library(dplyr)
library(xlsx)

metrics <- c('InDgr', 'OutDgr')  #, 'InDgr', 'OutDgr'

for (metric in metrics) {
  print(paste0("Starting analysis for metric: ", metric))  # Checkpoint 1
  
  outputFile <- paste0("c:/Users/j_m289/Pictures/phd/3. Data Analysis/studies/OKTOS/resting_te/analysis/brms/output/", tolower(metric), "_modelcomp.xlsx")

  ## ----Load RDS files -------------------------------------------------
  print("Loading model RDS files...")  # Checkpoint 2
  model1 <- readRDS(paste0("c:/Users/j_m289/Pictures/phd/3. Data Analysis/studies/OKTOS/resting_te/analysis/brms/models/", tolower(metric), "/", tolower(metric), "_model1.rds"))
  model2 <- readRDS(paste0("c:/Users/j_m289/Pictures/phd/3. Data Analysis/studies/OKTOS/resting_te/analysis/brms/models/", tolower(metric), "/", tolower(metric), "_model2.rds"))
  model3 <- readRDS(paste0("c:/Users/j_m289/Pictures/phd/3. Data Analysis/studies/OKTOS/resting_te/analysis/brms/models/", tolower(metric), "/", tolower(metric), "_model3.rds"))
  model4 <- readRDS(paste0("c:/Users/j_m289/Pictures/phd/3. Data Analysis/studies/OKTOS/resting_te/analysis/brms/models/", tolower(metric), "/", tolower(metric), "_model4.rds"))
  model5 <- readRDS(paste0("c:/Users/j_m289/Pictures/phd/3. Data Analysis/studies/OKTOS/resting_te/analysis/brms/models/", tolower(metric), "/", tolower(metric), "_model5.rds"))

  # ## ---- Assess fit (observed vs. predicted) -------------------------------------------------
  print("Calculating Bayes R2 for models...")  # Checkpoint 6
  
  # # Measure time for the first model
  start_time <- Sys.time()
  temp1 <- bayes_R2(model1) %>% as.data.frame()
  end_time <- Sys.time()
  time_taken_model1_bayesr2 <- end_time - start_time
  print(paste0("Time taken to calculate Bayes R2 for model1: ", time_taken_model1_bayesr2, ". Adjust expectations accordingly."))  # Time for model 1
  
  temp2 <- bayes_R2(model2) %>% as.data.frame()
  temp3 <- bayes_R2(model3) %>% as.data.frame()
  temp4 <- bayes_R2(model4) %>% as.data.frame()
  temp5 <- bayes_R2(model5) %>% as.data.frame()
  
  BayesR2 <- rbind(temp1, temp2, temp3, temp4, temp5)
  BayesR2 <- round(BayesR2, 4)
  
  row.names(BayesR2) <- paste("Model", 1:nrow(BayesR2))
  print("Bayes R2 calculations complete.")  # Checkpoint 7
  
  # # ---- Model comparison  -------------------------------------------------
  print("Performing model comparison...")  # Checkpoint 8
  loo <- loo(model1, model2, model3, model4, model5, compare = TRUE) # requires same numbers of observations
  modelcomparison <- loo$diffs %>% as.data.frame()
  print("Model comparison complete.")  # Checkpoint 9
  
  # ## ----- Data Export --------------
  print("Exporting results to Excel...")  # Checkpoint 10
  write.xlsx(BayesR2, outputFile, sheetName=paste0(metric,"_bayesR2"),  append=TRUE)
  write.xlsx(modelcomparison, outputFile, sheetName=paste0(metric,"_modelselection"),  append=TRUE)
  print("Results exported successfully.")  # Checkpoint 11

  gc()
  
  print(paste0("Analysis for metric ", metric, " complete."))  # Checkpoint 12
}
