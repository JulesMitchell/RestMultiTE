### MODEL SPECIFICATION AND VALIDATION ###
## ----pkgload, include = FALSE-------------------------------------------------
library(brms)
library(dplyr)
library(xlsx)

# Set output file name
outputFile <- "c:/Users/j_m289/Pictures/phd/3. Data Analysis/studies/OKTOS/resting_te/analysis/brms/output/asymmetry_finaloutput.xlsx"

## ----Load RDS files -------------------------------------------------
finalModel <- readRDS("c:/Users/j_m289/Pictures/phd/3. Data Analysis/studies/OKTOS/resting_te/analysis/brms/models/asymmetry/asymmetry_finalmodel.rds")

## ---- Add information criterion to model object  -------------------------------------------------
finalModel <- add_criterion(finalModel, c("waic", "loo"), save_psis = TRUE, overwrite = TRUE, file = 'c:/Users/j_m289/Pictures/phd/3. Data Analysis/studies/OKTOS/resting_te/analysis/brms/models/asymmetry/asymmetry_finalmodel')

## ---- Assess fit (observed vs. predicted) -------------------------------------------------
temp1 <- bayes_R2(finalModel) %>%as.data.frame()

loo <- finalModel$criteria$loo$estimates  %>% as.data.frame()
row.names(loo)[row.names(loo) == "model1"] <- finalModel$formula

loo(finalModel) # Copy pareto-k diagnostics

## ----- Data Export --------------
# Export tables
write.xlsx(temp1, outputFile, sheetName="BayesR2",  append=TRUE)
write.xlsx(loo, outputFile, sheetName="loo",  append=TRUE)