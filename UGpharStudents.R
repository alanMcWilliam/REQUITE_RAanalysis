###
### need patient demographics - adjust variables set
### need toxicity information
### need available pharmacy data
###


library(dplyr)
library(summarytools)

library(gtsummary)

# prostate data
prosFactor <- read.csv("C:/Users/alan_/Desktop/RAanalysis/REQUITEdata/Prostate/study/datasets/dataset5010.tsv", sep="\t", header=T)
View(prosFactor)

colnames(prosFactor)

prostatePharmacy <- prosFactor %>%
  select(SubjectId, ace_inhibitor, beta_blocker, other_antihypertensive_drug, on_statin, other_lipid_lowering_drugs, antidiabetic, pde5_inhibitor, sildenafil, X5alpha_reductase_inhibitor, alpha_blocker, anti_muscarinic, amiodarone, analgesics, antidepressant) %>%
  mutate(across(where(is.integer), as.factor))
summary(prostatePharmacy)
View(prostatePharmacy)

t <- summary(prostatePharmacy)
view(dfSummary(prostatePharmacy))

View(toxicityFiltered)

test <- merge(prostatePharmacy, toxicityFiltered, by = 'SubjectId')

tmp <- t(t)
knitr::kable(tmp)

test2 <- prostatePharmacy %>% 
  select(ace_inhibitor, beta_blocker, on_statin)
table <- tbl_summary(test2)
table

prostateComrob <- prosFactor %>%
  select(diabetes, ra, systemic_lupus_erythematosus, other_collagen_vascular_disease, hypertension, history_of_heart_disease, inflamatory_bowel_diventricular_disease, haemorrhoids, depression)%>%
  mutate(across(where(is.integer), as.factor))
summary(prostateComrob)

t <- summary(prostateComrob)
view(dfSummary(prostateComrob))

tmp <- t(t)
knitr::kable(tmp)

# lung data
lungFactor <- read.csv("C:/Users/alan_/Desktop/RAanalysis/REQUITEdata/Lung/study/datasets/dataset5022.tsv", sep="\t", header=T)
View(lungFactor)

colnames(lungFactor)

lungPharmacy <- lungFactor %>%
  select(hormone_replacement_therapy, on_statin, other_lipid_lowering_drugs, ace_inhibitor, other_antihypertensive_drug, amiodarone, antidiabetic, oral_steroids, analgesics, antidepressant, immunosuppressant) %>%
  mutate(across(where(is.integer), as.factor))
summary(lungPharmacy)

t <- summary(lungPharmacy)
view(dfSummary(lungPharmacy))

tmp <- t(t)
knitr::kable(tmp)

lungComrob <- lungFactor %>%
  select(diabetes, ra, systemic_lupus_erythematosus, other_collagen_vascular_disease, hypertension, history_of_heart_disease, depression,tuberculosis, copd) %>%
  mutate(across(where(is.integer), as.factor))
summary(lungComrob)

t <- summary(lungComrob)
view(dfSummary(lungComrob))

tmp <- t(t)
knitr::kable(tmp)

write.csv(test, "C:/Users/alan_/Desktop/REQUITE_polypharmacyData/test.csv")
