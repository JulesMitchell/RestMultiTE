### MODEL SPECIFICATION AND VALIDATION ###
library(brms)
library(dplyr)
library(xlsx)

outputFile <- "c:/Users/j_m289/Pictures/phd/3. Data Analysis/studies/OKTOS/resting_te/analysis/brms/output/geff_modelcomp.xlsx"

## ----Load RDS files -------------------------------------------------
model1 <- readRDS("c:/Users/j_m289/Pictures/phd/3. Data Analysis/studies/OKTOS/resting_te/analysis/brms/models/geff/geff_model1.rds")
model2 <- readRDS("c:/Users/j_m289/Pictures/phd/3. Data Analysis/studies/OKTOS/resting_te/analysis/brms/models/geff/geff_model2.rds")
model3 <- readRDS("c:/Users/j_m289/Pictures/phd/3. Data Analysis/studies/OKTOS/resting_te/analysis/brms/models/geff/geff_model3.rds")
model4 <- readRDS("c:/Users/j_m289/Pictures/phd/3. Data Analysis/studies/OKTOS/resting_te/analysis/brms/models/geff/geff_model4.rds")
model5 <- readRDS("c:/Users/j_m289/Pictures/phd/3. Data Analysis/studies/OKTOS/resting_te/analysis/brms/models/geff/geff_model5.rds")

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