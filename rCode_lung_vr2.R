

############################################################################################
####
#### A McWilliam
#### vr2 20th June 2022 updated for new lung analysis
#### 
#### late toxicity
####
############################################################################################


### gender, age, smoking status, concurrent chemotherapy, radiotherapy technique, FEV1, V20 Lung and V35 Esophagus
### smoking - never / Ex / Current
### chemo - yes / no
### RT technique - 3D-CRT / ARC / IMRT / Tomo / SBRT
### consider COPD

### cough, dyspnoea, pneumoitis, dysphagia, oesophagitis
### gender, age, smoking status, concurrent chemotherapy, radiotherapy technique, FEV1, V20 Lung and V35 Esophagus


### snps Ana Vega couldn't find in their datasets. Check impact
###rs187786174
###rs74984480
###rs2236668
###rs201408742
###rs5987194


############################################################################################
#### librarys needed
library(ggplot2)
library(dplyr)
library(summarytools)
library(lubridate)
library(ggpubr)
library(stringr)


############################################################################################
### Data cleaning and calculating the STAT


### read in prostate toxicities
### also need to know data of radiotherapy
Tox <- read.csv("C:/Users/alan_/Desktop/RAanalysis/REQUITEdata/Lung/study/datasets/dataset5025.tsv", sep="\t", header=T )
View(Tox)
Treat <- read.csv("C:/Users/alan_/Desktop/RAanalysis/REQUITEdata/Lung/study/datasets/dataset5029.tsv", sep="\t", header=T)
View(Treat)
Factor <- read.csv("C:/Users/alan_/Desktop/RAanalysis/REQUITEdata/Lung/study/datasets/dataset5022.tsv", sep="\t", header=T)
View(Factor)

lungPRO <- read.csv("C:/Users/alan_/Desktop/RAanalysis/REQUITEdata/Lung/study/datasets/dataset5004.tsv", sep="\t", header=T)
View(lungPRO)

### select ID and radiotherapy start data
radiotherapyStart <- Treat %>%
  select(SubjectId, radiotherapy_start_date)


### select out toxicity data of interest for STAT
#lungPRO$L5a_q05_cough, lungPRO$L5a_q06_cough_blood, lungPRO$L5a_q04_short_breath, lungPRO$L5a_q02_problems_swallowing, lungPRO$L5a_q07_lost_appetite ,lungPRO$L5a_q03_pain_chest
lungTox <- Tox %>%
  select(SubjectId, date, cough, dyspnoea, pneumonitis, dysphagia, esophagitis)


toxicity <- merge(lungTox, radiotherapyStart, by = 'SubjectId')
toxicity$monthStartTreat <- interval(ymd(as.Date(toxicity$radiotherapy_start_date)), ymd(as.Date(toxicity$date)))
toxicity$monthStartTreat = toxicity$monthStartTreat %/% months(1)

View(toxicity)

### 

### need to set a time point for selecting highest toxicity score and filter
### keep highest recorded score for each toxicity 
months = 500 ## acute toxicity - 90 days acute set months = 4 minMonths = 1  late 500 (big number) min months 3
minMonths = 3

## identify patients with baseline values
baselineCounts <- toxicity %>%
  group_by(SubjectId) %>%
  filter(monthStartTreat == 0) %>%
  select(SubjectId)


## filter by baseline and subtract first entry
toxicitySubtract <- toxicity %>%
  select(SubjectId, monthStartTreat, cough, dyspnoea, pneumonitis, dysphagia, esophagitis) %>%
  filter(SubjectId %in% baselineCounts$SubjectId) %>%
  group_by(SubjectId) %>%
  mutate(cough_diff = cough - first(cough), 
         dyspnoea_diff = dyspnoea - first(dyspnoea), 
         pneumonitis_diff = pneumonitis - first(pneumonitis), 
         dysphagia_diff = dysphagia - first(dysphagia), 
         esophagitis_diff = esophagitis - first(esophagitis)) %>%
  filter(monthStartTreat > minMonths) %>%
  filter(monthStartTreat < months) %>%
  summarise(maxMonths = max(monthStartTreat), 
            maxCough = max(cough), 
            maxDyspnoea = max(dyspnoea), 
            maxPneumonitis = max(pneumonitis), 
            maxDysphagia = max(dysphagia), 
            maxEsophagitis = max(esophagitis)) %>%
  select(SubjectId, maxMonths, maxCough, maxDyspnoea, maxPneumonitis, maxDysphagia, maxEsophagitis)

toxicitySubtract[toxicitySubtract < 0] <- 0 
View(toxicitySubtract)


## patients without baseline
toxicityNoBaseline <- toxicity %>%
  filter(!SubjectId %in% baselineCounts$SubjectId) %>%
  filter(monthStartTreat > minMonths) %>%
  filter(monthStartTreat < months) %>%
  group_by(SubjectId) %>%
  summarise(maxMonths = max(monthStartTreat), 
            maxCough = max(cough), 
            maxDyspnoea = max(dyspnoea), 
            maxPneumonitis = max(pneumonitis), 
            maxDysphagia = max(dysphagia), 
            maxEsophagitis = max(esophagitis)) %>%
  select(SubjectId, maxMonths, maxCough, maxDyspnoea, maxPneumonitis, maxDysphagia, maxEsophagitis)
View(toxicityNoBaseline)

## join both together again
toxicityFiltered <- rbind(toxicitySubtract, toxicityNoBaseline)
View(toxicityFiltered)


###create temp data frame to display summary stats
toxicity_summaryStats <- toxicityFiltered %>%
  select( maxCough, maxDyspnoea, maxPneumonitis, maxDysphagia, maxEsophagitis)

stview(dfSummary(toxicity_summaryStats))


t <- toxicityFiltered[ ,3:ncol(toxicityFiltered)]
col_Mean <- t %>% summarise_if(is.numeric, mean, na.rm = TRUE)
col_SD <- t %>% summarise_if(is.numeric, sd, na.rm = TRUE)
STAT <- matrix(NA, nrow = nrow(t), ncol = 1)

for(i in 1:nrow(t)){
  tmp = 0
  count = 0
  for(j in 1:length(t)){
    if(!is.na(t[i,j])){
      tmp = tmp + (t[i,j] - col_Mean[j])/col_SD[j]
      count = count + 1
    }
  }
  tmp = tmp/count
  STAT[i] = tmp
}

View(STAT)

### join STAT back on to data 
toxicityFilteredSTAT <- bind_cols(toxicityFiltered, t(as.data.frame(STAT)))
names(toxicityFilteredSTAT)[ncol(toxicityFilteredSTAT)] <- "STAT"

toxicityFilteredSTAT <- toxicityFilteredSTAT[!is.nan(toxicityFilteredSTAT$STAT), ]
View(toxicityFilteredSTAT)

### select other clinical variables for inclusion in analysis


### save csv of data wkith STAT ready for analysis.
write.csv(toxicityFilteredSTAT, "C:/Users/alan_/Desktop/RAanalysis/REQUITEdata/processed/lung_acuteSTAT.csv")
#write.csv(toxicityFilteredSTAT, "C:/Users/alan_/Desktop/RAanalysis/REQUITEdata/processed/lung_lateSTAT.csv")

############################################################################################
### summary of patient characteristics
### summary STAT and plot histogram
###

summary(toxicityFilteredSTAT$STAT)

ggplot(data=toxicityFilteredSTAT, aes(STAT)) + 
  geom_histogram(breaks=seq(-1,3, by = 0.2),
                 col = "skyblue", fill = "lightblue") +
  labs(title = i, x = "STAT" ) +
  theme(panel.background = element_blank())


#ggplot(data=toxicityFilteredSTAT, aes(log10(STAT))) + 
#  geom_histogram(breaks=seq(-3,3, by = 0.1),
#                 col = "skyblue", fill = "lightblue") +
#  labs(title = i, x = "STAT" ) +
#  theme(panel.background = element_blank())


### save plots - change name
#ggsave("C:/Users/alan_/Desktop/RAanalysis/REQUITEdata/processed/figures/STATprostate.jpg")


##########################################################################################
### load in PRS + wPRS

alanPRS <- read.csv("C:/Users/alan_/Desktop/RAanalysis/calcPRS/PRS_all_NEW.csv", header = F)
alanPRS <- t(alanPRS)
alanPRS <- alanPRS[-1,]

colnames(alanPRS) <- c("SampleID", "prs_alan")
alanPRS <- as.data.frame(alanPRS)
alanPRS$SampleID <- as.numeric(alanPRS$SampleID)
alanPRS$prs_alan <- as.numeric(alanPRS$prs_alan)
#View(alanPRS)

alanWPRS <- read.csv("C:/Users/alan_/Desktop/RAanalysis/calcPRS/wPRS_NEW.csv", header = F)
alanWPRS <- t(alanWPRS)
alanWPRS <- alanWPRS[-1,]

colnames(alanWPRS) <- c("SampleID", "wprs_alan")
alanWPRS <- as.data.frame(alanWPRS)
alanWPRS$SampleID <- as.numeric(alanWPRS$SampleID)
alanWPRS$wprs_alan <- as.numeric(alanWPRS$wprs_alan)
#View(alanWPRS)

PRS_all <- merge(alanPRS, alanWPRS, by = "SampleID")
View(PRS_all)

ggplot(data = PRS_all) +
  geom_histogram(aes(x = prs_alan), 
                 binwidth = 1, col = "skyblue", fill = "lightblue") +
  labs(title = "", x = "prs") +
  theme(panel.background = element_blank())


ggplot(data = PRS_all) +
  geom_histogram(aes(x = wprs_alan), 
                 binwidth = 0.1, col = "skyblue", fill = "lightblue") +
  labs(title = "", x = "wprs") +
  theme(panel.background = element_blank())

summary(PRS_all$prs_alan)
summary(PRS_all$wprs_alan)

## need to link back to patientID
sampleIDlink <- read.csv("C:/Users/alan_/Desktop/RAanalysis/REQUITEdata/Lung/study/datasets/dataset5039.tsv", sep="\t", header=T)

##t <- merge(sampleIDlink, toxicityCROFilteredSTAT, by = "SubjectId")
STAT_prs <- merge(PRS_all, sampleIDlink, by = "SampleID")
STAT_prs <- merge(STAT_prs, toxicityFilteredSTAT, by = "SubjectId")
View(STAT_prs)


#t <- glm(STAT~prs_precentile_alan, data = STAT_prs)
#STAT_prs$prs_precentile_alan <- STAT_prs$prs_alan > quantile(STAT_prs$prs_alan, c(.95)) 
#summary(t)

##########################################################################################
### need to select other patient factors
### adjusting for androgen-deprivation therapy, prior prostatectomy, age at treatment, and total BED
### + diabetes? 
### BED for acute tox - 3 also try 10 for normal tissue???

patFEV <- Tox %>%
  select(SubjectId, date, fev1_litres) %>%
  na.omit() %>%
  group_by(SubjectId) %>%
  mutate(fev = first(fev1_litres)) %>%
  filter(row_number()==1) %>%
  filter(fev <100) %>%
  select(SubjectId, fev)
summary(patFEV)

patFactors <- Factor %>%
  select(SubjectId, gender, age_at_radiotherapy_start_yrs, smoker, diabetes, ra, copd)

patTreat <- Treat %>%
  select( SubjectId, radiotherapy_number_of_fractions, radiotherapy_total_dose_Gy, chemotherapy, radiotherapy_technique, radiotherapy_v20_lung_pc, radiotherapy_v35_oesophagus_pc) ##, p3radical_prostatectomy, p3hormone_therapy, p3hormone_therapy_length_months)


### need to calculate the BED prescribed
### BED = D x (1 + [d / (??/??)])
alphaBeta <- 10 #3
patTreat$doseFraction <- patTreat$radiotherapy_total_dose_Gy/patTreat$radiotherapy_number_of_fractions
patTreat$doseBED <-   patTreat$radiotherapy_total_dose_Gy  * (1 + (patTreat$doseFraction/alphaBeta))



STAT_prs_factors <- merge(STAT_prs, patFactors, by = "SubjectId")
STAT_prs_factors <- merge(STAT_prs_factors, patTreat, by = "SubjectId")
STAT_prs_factors <- merge(STAT_prs_factors, patFEV, by = "SubjectId")


### clean the data
#STAT_prs_factors <- STAT_prs_factors %>%
#  filter(CancerType == 'Prostate')

View(STAT_prs_factors)

view(dfSummary(STAT_prs_factors))



patChra <- merge(patFactors, patTreat, by = "SubjectId")
patChra <- merge(patChra, patFEV, by = "SubjectId")
view(dfSummary(patChra))
  
##########################################################################################
### analysing PRS and wPRS
summary(STAT_prs_factors$STAT)

STAT_prs_factors <- STAT_prs_factors %>% mutate_all(~ifelse(is.nan(.), NA, .))


### lost some code here - check correct
t <- glm(STAT~prs_alan, data = STAT_prs_factors)
summary(t)
confint(t)

t <- glm(STAT~wprs_alan, data = STAT_prs_factors)
summary(t)
confint(t)

STAT_prs_factors$prs_precentile_alan <- STAT_prs_factors$prs_alan > quantile(STAT_prs_factors$prs_alan, c(.90)) 
t <- glm(STAT~prs_precentile_alan, data = STAT_prs_factors)
summary(t)
confint(t)

STAT_prs_factors$wprs_precentile_alan <- STAT_prs_factors$wprs_alan > quantile(STAT_prs_factors$wprs_alan, c(.90)) 
t <- glm(STAT~wprs_precentile_alan, data = STAT_prs_factors)
summary(t)
confint(t)

t <- glm(STAT~prs_alan + factor(gender) + age_at_radiotherapy_start_yrs + factor(smoker) + factor(chemotherapy) + fev + radiotherapy_v20_lung_pc + radiotherapy_v35_oesophagus_pc + doseBED + factor(copd), data = STAT_prs_factors)
summary(t) 
confint(t)
t <- glm(STAT~wprs_alan + factor(gender) + age_at_radiotherapy_start_yrs + factor(smoker) + factor(chemotherapy) + fev + radiotherapy_v20_lung_pc + radiotherapy_v35_oesophagus_pc + doseBED + factor(copd), data = STAT_prs_factors)
summary(t) 
confint(t)

t <- glm(STAT~prs_precentile_alan + factor(gender) + age_at_radiotherapy_start_yrs + factor(smoker) + factor(chemotherapy) + factor(radiotherapy_technique) + fev + radiotherapy_v20_lung_pc + radiotherapy_v35_oesophagus_pc, data = STAT_prs_factors)
summary(t) 
confint(t)
t <- glm(STAT~wprs_precentile_alan + factor(gender) + age_at_radiotherapy_start_yrs + factor(smoker) + factor(chemotherapy) + fev + radiotherapy_v20_lung_pc + radiotherapy_v35_oesophagus_pc + doseBED + factor(copd), data = STAT_prs_factors)
summary(t) 
confint(t)

t <- glm(STAT~prs_precentile_alan + factor(gender) + age_at_radiotherapy_start_yrs + factor(smoker) + factor(chemotherapy) + fev + radiotherapy_v20_lung_pc + radiotherapy_v35_oesophagus_pc + doseBED + factor(copd), data = STAT_prs_factors)
summary(t) 




t <- glm(STAT~prs_precentile_alan + factor(gender) + age_at_radiotherapy_start_yrs + factor(smoker) + factor(chemotherapy) + fev + radiotherapy_v20_lung_pc + radiotherapy_v35_oesophagus_pc + doseBED + factor(copd), data = STAT_prs_factors)
summary(t) 

##### individual endpoints

## maxCough
#prs
tt1 <- glm(maxCough~prs_alan + factor(gender) + age_at_radiotherapy_start_yrs + factor(smoker) + factor(chemotherapy) + fev + radiotherapy_v20_lung_pc + radiotherapy_v35_oesophagus_pc + doseBED + factor(copd), data = STAT_prs_factors)
summary(tt1) 
ttt1 <- confint(tt1, level = 0.95)
a1 <- paste(format(round(ttt1[2,1], 3), nsmall = 3), format(round(ttt1[2,2], 3), nsmall = 3), sep=', ')
aa1<- paste('(',a1,')', sep = '')
aaa1 <- paste(format(round(as.numeric(summary(tt1)$coeff[2]), 3), nsmall = 3), aa1, format(round(summary(tt1)$coeff[2,4], 2), nsmall = 2), sep = ' ')


#wprs
tt1 <- glm(maxCough~wprs_alan + factor(gender) + age_at_radiotherapy_start_yrs + factor(smoker) + factor(chemotherapy) + fev + radiotherapy_v20_lung_pc + radiotherapy_v35_oesophagus_pc + doseBED + factor(copd), data = STAT_prs_factors)
summary(tt1) 
ttt1 <- confint(tt1, level = 0.95)
a2 <- paste(format(round(ttt1[2,1], 3), nsmall = 3), format(round(ttt1[2,2], 3), nsmall = 3), sep=', ')
aa2<- paste('(',a2,')', sep = '')
aaa2 <- paste(format(round(as.numeric(summary(tt1)$coeff[2]), 3), nsmall = 3), aa2, format(round(summary(tt1)$coeff[2,4], 2), nsmall = 2), sep = ' ')

aaa1
aaa2

#prs percentile
tt1 <- glm(maxCough~prs_precentile_alan + factor(gender) + age_at_radiotherapy_start_yrs + factor(smoker) + factor(chemotherapy) + fev + radiotherapy_v20_lung_pc + radiotherapy_v35_oesophagus_pc + doseBED + factor(copd), data = STAT_prs_factors)
summary(tt1) 
ttt1 <- confint(tt1, level = 0.95)
a1 <- paste(format(round(ttt1[2,1], 3), nsmall = 3), format(round(ttt1[2,2], 3), nsmall = 3), sep=', ')
aa1<- paste('(',a1,')', sep = '')
aaa1 <- paste(format(round(as.numeric(summary(tt1)$coeff[2]), 3), nsmall = 3), aa1, format(round(summary(tt1)$coeff[2,4], 2), nsmall = 2), sep = ' ')

#wprs percentile
tt1 <- glm(maxCough~wprs_precentile_alan + factor(gender) + age_at_radiotherapy_start_yrs + factor(smoker) + factor(chemotherapy) + fev + radiotherapy_v20_lung_pc + radiotherapy_v35_oesophagus_pc + doseBED + factor(copd), data = STAT_prs_factors)
summary(tt1) 
ttt1 <- confint(tt1, level = 0.95)
a2 <- paste(format(round(ttt1[2,1], 3), nsmall = 3), format(round(ttt1[2,2], 3), nsmall = 3), sep=', ')
aa2<- paste('(',a2,')', sep = '')
aaa2 <- paste(format(round(as.numeric(summary(tt1)$coeff[2]), 3), nsmall = 3), aa2, format(round(summary(tt1)$coeff[2,4], 2), nsmall = 2), sep = ' ')

aaa1
aaa2

## maxDyspnoea
#prs
tt1 <- glm(maxDyspnoea~prs_alan + factor(gender) + age_at_radiotherapy_start_yrs + factor(smoker) + factor(chemotherapy) + fev + radiotherapy_v20_lung_pc + radiotherapy_v35_oesophagus_pc + doseBED + factor(copd), data = STAT_prs_factors)
summary(tt1) 
ttt1 <- confint(tt1, level = 0.95)
a1 <- paste(format(round(ttt1[2,1], 3), nsmall = 3), format(round(ttt1[2,2], 3), nsmall = 3), sep=', ')
aa1<- paste('(',a1,')', sep = '')
aaa1 <- paste(format(round(as.numeric(summary(tt1)$coeff[2]), 3), nsmall = 3), aa1, format(round(summary(tt1)$coeff[2,4], 2), nsmall = 2), sep = ' ')


#wprs
tt1 <- glm(maxDyspnoea~wprs_alan + factor(gender) + age_at_radiotherapy_start_yrs + factor(smoker) + factor(chemotherapy) + fev + radiotherapy_v20_lung_pc + radiotherapy_v35_oesophagus_pc + doseBED + factor(copd), data = STAT_prs_factors)
summary(tt1) 
ttt1 <- confint(tt1, level = 0.95)
a2 <- paste(format(round(ttt1[2,1], 3), nsmall = 3), format(round(ttt1[2,2], 3), nsmall = 3), sep=', ')
aa2<- paste('(',a2,')', sep = '')
aaa2 <- paste(format(round(as.numeric(summary(tt1)$coeff[2]), 3), nsmall = 3), aa2, format(round(summary(tt1)$coeff[2,4], 2), nsmall = 2), sep = ' ')

aaa1
aaa2

#prs percentile
tt1 <- glm(maxDyspnoea~prs_precentile_alan + factor(gender) + age_at_radiotherapy_start_yrs + factor(smoker) + factor(chemotherapy) + fev + radiotherapy_v20_lung_pc + radiotherapy_v35_oesophagus_pc + doseBED + factor(copd), data = STAT_prs_factors)
summary(tt1) 
ttt1 <- confint(tt1, level = 0.95)
a1 <- paste(format(round(ttt1[2,1], 3), nsmall = 3), format(round(ttt1[2,2], 3), nsmall = 3), sep=', ')
aa1<- paste('(',a1,')', sep = '')
aaa1 <- paste(format(round(as.numeric(summary(tt1)$coeff[2]), 3), nsmall = 3), aa1, format(round(summary(tt1)$coeff[2,4], 2), nsmall = 2), sep = ' ')


#wprs percentile
tt1 <- glm(maxDyspnoea~wprs_precentile_alan + factor(gender) + age_at_radiotherapy_start_yrs + factor(smoker) + factor(chemotherapy) + fev + radiotherapy_v20_lung_pc + radiotherapy_v35_oesophagus_pc + doseBED + factor(copd), data = STAT_prs_factors)
summary(tt1) 
ttt1 <- confint(tt1, level = 0.95)
a2 <- paste(format(round(ttt1[2,1], 3), nsmall = 3), format(round(ttt1[2,2], 3), nsmall = 3), sep=', ')
aa2<- paste('(',a2,')', sep = '')
aaa2 <- paste(format(round(as.numeric(summary(tt1)$coeff[2]), 3), nsmall = 3), aa2, format(round(summary(tt1)$coeff[2,4], 2), nsmall = 2), sep = ' ')

aaa1
aaa2

##maxPneumonitis
#prs
tt1 <- glm(maxPneumonitis~prs_alan + factor(gender) + age_at_radiotherapy_start_yrs + factor(smoker) + factor(chemotherapy) + fev + radiotherapy_v20_lung_pc + radiotherapy_v35_oesophagus_pc + doseBED + factor(copd), data = STAT_prs_factors)
summary(tt1) 
ttt1 <- confint(tt1, level = 0.95)
a1 <- paste(format(round(ttt1[2,1], 3), nsmall = 3), format(round(ttt1[2,2], 3), nsmall = 3), sep=', ')
aa1<- paste('(',a1,')', sep = '')
aaa1 <- paste(format(round(as.numeric(summary(tt1)$coeff[2]), 3), nsmall = 3), aa1, format(round(summary(tt1)$coeff[2,4], 2), nsmall = 2), sep = ' ')


#wprs
tt1 <- glm(maxPneumonitis~wprs_alan + factor(gender) + age_at_radiotherapy_start_yrs + factor(smoker) + factor(chemotherapy) + fev + radiotherapy_v20_lung_pc + radiotherapy_v35_oesophagus_pc + doseBED + factor(copd), data = STAT_prs_factors)
summary(tt1) 
ttt1 <- confint(tt1, level = 0.95)
a2 <- paste(format(round(ttt1[2,1], 3), nsmall = 3), format(round(ttt1[2,2], 3), nsmall = 3), sep=', ')
aa2<- paste('(',a2,')', sep = '')
aaa2 <- paste(format(round(as.numeric(summary(tt1)$coeff[2]), 3), nsmall = 3), aa2, format(round(summary(tt1)$coeff[2,4], 2), nsmall = 2), sep = ' ')

aaa1
aaa2

#prs percentile
tt1 <- glm(maxPneumonitis~prs_precentile_alan + factor(gender) + age_at_radiotherapy_start_yrs + factor(smoker) + factor(chemotherapy) + fev + radiotherapy_v20_lung_pc + radiotherapy_v35_oesophagus_pc + doseBED + factor(copd), data = STAT_prs_factors)
summary(tt1) 
ttt1 <- confint(tt1, level = 0.95)
a1 <- paste(format(round(ttt1[2,1], 3), nsmall = 3), format(round(ttt1[2,2], 3), nsmall = 3), sep=', ')
aa1<- paste('(',a1,')', sep = '')
aaa1 <- paste(format(round(as.numeric(summary(tt1)$coeff[2]), 3), nsmall = 3), aa1, format(round(summary(tt1)$coeff[2,4], 2), nsmall = 2), sep = ' ')


#wprs percentile
tt1 <- glm(maxPneumonitis~wprs_precentile_alan + factor(gender) + age_at_radiotherapy_start_yrs + factor(smoker) + factor(chemotherapy) + fev + radiotherapy_v20_lung_pc + radiotherapy_v35_oesophagus_pc + doseBED + factor(copd), data = STAT_prs_factors)
summary(tt1) 
ttt1 <- confint(tt1, level = 0.95)
a2 <- paste(format(round(ttt1[2,1], 3), nsmall = 3), format(round(ttt1[2,2], 3), nsmall = 3), sep=', ')
aa2<- paste('(',a2,')', sep = '')
aaa2 <- paste(format(round(as.numeric(summary(tt1)$coeff[2]), 3), nsmall = 3), aa2, format(round(summary(tt1)$coeff[2,4], 2), nsmall = 2), sep = ' ')

aaa1
aaa2

##maxDysphagia
#prs
tt1 <- glm(maxDysphagia~prs_alan + factor(gender) + age_at_radiotherapy_start_yrs + factor(smoker) + factor(chemotherapy) + fev + radiotherapy_v20_lung_pc + radiotherapy_v35_oesophagus_pc + doseBED + factor(copd), data = STAT_prs_factors)
summary(tt1) 
ttt1 <- confint(tt1, level = 0.95)
a1 <- paste(format(round(ttt1[2,1], 3), nsmall = 3), format(round(ttt1[2,2], 3), nsmall = 3), sep=', ')
aa1<- paste('(',a1,')', sep = '')
aaa1 <- paste(format(round(as.numeric(summary(tt1)$coeff[2]), 3), nsmall = 3), aa1, format(round(summary(tt1)$coeff[2,4], 2), nsmall = 2), sep = ' ')


#wprs
tt1 <- glm(maxDysphagia~wprs_alan + factor(gender) + age_at_radiotherapy_start_yrs + factor(smoker) + factor(chemotherapy) + fev + radiotherapy_v20_lung_pc + radiotherapy_v35_oesophagus_pc + doseBED + factor(copd), data = STAT_prs_factors)
summary(tt1) 
ttt1 <- confint(tt1, level = 0.95)
a2 <- paste(format(round(ttt1[2,1], 3), nsmall = 3), format(round(ttt1[2,2], 3), nsmall = 3), sep=', ')
aa2<- paste('(',a2,')', sep = '')
aaa2 <- paste(format(round(as.numeric(summary(tt1)$coeff[2]), 3), nsmall = 3), aa2, format(round(summary(tt1)$coeff[2,4], 2), nsmall = 2), sep = ' ')

aaa1
aaa2

#prs percentile
tt1 <- glm(maxDysphagia~prs_precentile_alan + factor(gender) + age_at_radiotherapy_start_yrs + factor(smoker) + factor(chemotherapy) + fev + radiotherapy_v20_lung_pc + radiotherapy_v35_oesophagus_pc + doseBED + factor(copd), data = STAT_prs_factors)
summary(tt1) 
ttt1 <- confint(tt1, level = 0.95)
a1 <- paste(format(round(ttt1[2,1], 3), nsmall = 3), format(round(ttt1[2,2], 3), nsmall = 3), sep=', ')
aa1<- paste('(',a1,')', sep = '')
aaa1 <- paste(format(round(as.numeric(summary(tt1)$coeff[2]), 3), nsmall = 3), aa1, format(round(summary(tt1)$coeff[2,4], 2), nsmall = 2), sep = ' ')


#wprs percentile
tt1 <- glm(maxDysphagia~wprs_precentile_alan + factor(gender) + age_at_radiotherapy_start_yrs + factor(smoker) + factor(chemotherapy) + fev + radiotherapy_v20_lung_pc + radiotherapy_v35_oesophagus_pc + doseBED + factor(copd), data = STAT_prs_factors)
summary(tt1) 
ttt1 <- confint(tt1, level = 0.95)
a2 <- paste(format(round(ttt1[2,1], 3), nsmall = 3), format(round(ttt1[2,2], 3), nsmall = 3), sep=', ')
aa2<- paste('(',a2,')', sep = '')
aaa2 <- paste(format(round(as.numeric(summary(tt1)$coeff[2]), 3), nsmall = 3), aa2, format(round(summary(tt1)$coeff[2,4], 2), nsmall = 2), sep = ' ')

aaa1
aaa2


## maxEsophagitis

#prs
tt1 <- glm(maxEsophagitis~prs_alan + factor(gender) + age_at_radiotherapy_start_yrs + factor(smoker) + factor(chemotherapy) + fev + radiotherapy_v20_lung_pc + radiotherapy_v35_oesophagus_pc + doseBED + factor(copd), data = STAT_prs_factors)
summary(tt1) 
ttt1 <- confint(tt1, level = 0.95)
a1 <- paste(format(round(ttt1[2,1], 3), nsmall = 3), format(round(ttt1[2,2], 3), nsmall = 3), sep=', ')
aa1<- paste('(',a1,')', sep = '')
aaa1 <- paste(format(round(as.numeric(summary(tt1)$coeff[2]), 3), nsmall = 3), aa1, format(round(summary(tt1)$coeff[2,4], 2), nsmall = 2), sep = ' ')


#wprs
tt1 <- glm(maxEsophagitis~wprs_alan + factor(gender) + age_at_radiotherapy_start_yrs + factor(smoker) + factor(chemotherapy) + fev + radiotherapy_v20_lung_pc + radiotherapy_v35_oesophagus_pc + doseBED + factor(copd), data = STAT_prs_factors)
summary(tt1) 
ttt1 <- confint(tt1, level = 0.95)
a2 <- paste(format(round(ttt1[2,1], 3), nsmall = 3), format(round(ttt1[2,2], 3), nsmall = 3), sep=', ')
aa2<- paste('(',a2,')', sep = '')
aaa2 <- paste(format(round(as.numeric(summary(tt1)$coeff[2]), 3), nsmall = 3), aa2, format(round(summary(tt1)$coeff[2,4], 2), nsmall = 2), sep = ' ')

aaa1
aaa2

#prs percentile
tt1 <- glm(maxEsophagitis~prs_precentile_alan + factor(gender) + age_at_radiotherapy_start_yrs + factor(smoker) + factor(chemotherapy) + fev + radiotherapy_v20_lung_pc + radiotherapy_v35_oesophagus_pc + doseBED + factor(copd), data = STAT_prs_factors)
summary(tt1) 
ttt1 <- confint(tt1, level = 0.95)
a1 <- paste(format(round(ttt1[2,1], 3), nsmall = 3), format(round(ttt1[2,2], 3), nsmall = 3), sep=', ')
aa1<- paste('(',a1,')', sep = '')
aaa1 <- paste(format(round(as.numeric(summary(tt1)$coeff[2]), 3), nsmall = 3), aa1, format(round(summary(tt1)$coeff[2,4], 2), nsmall = 2), sep = ' ')


#wprs percentile
tt1 <- glm(maxEsophagitis~wprs_precentile_alan + factor(gender) + age_at_radiotherapy_start_yrs + factor(smoker) + factor(chemotherapy) + fev + radiotherapy_v20_lung_pc + radiotherapy_v35_oesophagus_pc + doseBED + factor(copd), data = STAT_prs_factors)
summary(tt1) 
ttt1 <- confint(tt1, level = 0.95)
a2 <- paste(format(round(ttt1[2,1], 3), nsmall = 3), format(round(ttt1[2,2], 3), nsmall = 3), sep=', ')
aa2<- paste('(',a2,')', sep = '')
aaa2 <- paste(format(round(as.numeric(summary(tt1)$coeff[2]), 3), nsmall = 3), aa2, format(round(summary(tt1)$coeff[2,4], 2), nsmall = 2), sep = ' ')

aaa1
aaa2








##########################################

t <- glm(maxCough~prs_precentile_alan + factor(gender) + age_at_radiotherapy_start_yrs + factor(smoker) + factor(chemotherapy) + fev + radiotherapy_v20_lung_pc + radiotherapy_v35_oesophagus_pc + doseBED + factor(copd), data = STAT_prs_factors)
summary(t) 
t <- glm(maxDyspnoea~prs_precentile_alan + factor(gender) + age_at_radiotherapy_start_yrs + factor(smoker) + factor(chemotherapy) + fev + radiotherapy_v20_lung_pc + radiotherapy_v35_oesophagus_pc + doseBED + factor(copd), data = STAT_prs_factors)
summary(t) 
t <- glm(maxPneumonitis~prs_precentile_alan + factor(gender) + age_at_radiotherapy_start_yrs + factor(smoker) + factor(chemotherapy) + fev + doseBED + radiotherapy_v20_lung_pc + radiotherapy_v35_oesophagus_pc + factor(copd), data = STAT_prs_factors)
summary(t) 
t <- glm(maxDysphagia~prs_precentile_alan + factor(gender) + age_at_radiotherapy_start_yrs + factor(smoker) + factor(chemotherapy) + fev + radiotherapy_v20_lung_pc + radiotherapy_v35_oesophagus_pc + doseBED + factor(copd), data = STAT_prs_factors)
summary(t) 
t <- glm(maxEsophagitis~prs_precentile_alan + factor(gender) + age_at_radiotherapy_start_yrs + factor(smoker) + factor(chemotherapy) + fev + radiotherapy_v20_lung_pc + radiotherapy_v35_oesophagus_pc + doseBED + factor(copd), data = STAT_prs_factors)
summary(t) 



t <- glm(maxCough~wprs_precentile_alan + factor(gender) + age_at_radiotherapy_start_yrs + factor(smoker) + factor(chemotherapy) + fev + radiotherapy_v20_lung_pc + radiotherapy_v35_oesophagus_pc + doseBED + factor(copd), data = STAT_prs_factors)
summary(t) 
t <- glm(maxDyspnoea~wprs_precentile_alan + factor(gender) + age_at_radiotherapy_start_yrs + factor(smoker) + factor(chemotherapy) + fev + radiotherapy_v20_lung_pc + radiotherapy_v35_oesophagus_pc + doseBED + factor(copd), data = STAT_prs_factors)
summary(t) 
t <- glm(maxPneumonitis~wprs_precentile_alan + factor(gender) + age_at_radiotherapy_start_yrs + factor(smoker) + factor(chemotherapy) + fev + radiotherapy_v20_lung_pc + radiotherapy_v35_oesophagus_pc + doseBED + factor(copd), data = STAT_prs_factors)
summary(t) 
t <- glm(maxDysphagia~wprs_precentile_alan + factor(gender) + age_at_radiotherapy_start_yrs + factor(smoker) + factor(chemotherapy) + fev + radiotherapy_v20_lung_pc + radiotherapy_v35_oesophagus_pc + doseBED + factor(copd), data = STAT_prs_factors)
summary(t) 
t <- glm(maxEsophagitis~wprs_precentile_alan + factor(gender) + age_at_radiotherapy_start_yrs + factor(smoker) + factor(chemotherapy) + fev + radiotherapy_v20_lung_pc + radiotherapy_v35_oesophagus_pc + doseBED + factor(copd), data = STAT_prs_factors)
summary(t) 






##########################################################################################
##########################################################################################
##########################################################################################


#######################################################
#######################################################
### look at patients with Rhematoid arthritis
summary(factor(STAT_prs_factors$ra))
tapply(STAT_prs_factors$wprs_alan, STAT_prs_factors$ra, summary)
tapply(STAT_prs_factors$prs_alan, STAT_prs_factors$ra, summary)

STAT_prs_factors$ra <- as.factor(STAT_prs_factors$ra)
ggplot(STAT_prs_factors, aes(x = ra, y = prs_alan)) + 
  geom_boxplot() + 
  theme_classic() +
  stat_compare_means(method = "wilcox.test", aes(group = ra, label = paste0("p = ",..p.format..)), label.x = 1.4, label.y = 6, size = 6)

ggplot(STAT_prs_factors, aes(x = ra, y = wprs_alan)) + 
  geom_boxplot() + 
  theme_classic() +
  stat_compare_means(method = "wilcox.test", aes(group = ra, label = paste0("p = ",..p.format..)), label.x = 1.4, label.y = 6, size = 6)





##########################################################################################
##### work with residuals

t2 <- glm(STAT~ factor(gender) + age_at_radiotherapy_start_yrs + factor(smoker) + factor(chemotherapy) + fev + radiotherapy_v20_lung_pc + radiotherapy_v35_oesophagus_pc + doseBED + factor(copd) + diabetes, data = STAT_prs_factors)
summary(t2) 
AIC(t2)

resid(t2)
plot(density(resid(t2))) #A density plot
qqnorm(resid(t2)) # A quantile normal plot - good for checking normality
qqline(resid(t2))


residualsTmp <- t2$residuals 
residuals_present <- intersect(names(residualsTmp),rownames(STAT_prs_factors))
STAT_residuals <- STAT_prs_factors[residuals_present,]
STAT_residuals$residual<-residualsTmp

View(STAT_residuals)

t_residual <- glm(wprs_alan~residual, data = STAT_residuals)
summary(t_residual)
t_residual <- glm(prs_alan~residual, data = STAT_residuals)
summary(t_residual)

t_residual <- glm(prs_precentile_alan~residual, data = STAT_residuals)
summary(t_residual)
t_residual <- glm(wprs_precentile_alan~residual, data = STAT_residuals)
summary(t_residual)


### read in genedoses and format dataframe
genedose <- read.csv("C:/Users/alan_/Desktop/RAanalysis/calcPRS/genedoses.csv", header = F)
genedose <- t(genedose)
genedose <- genedose[-1,]
genedose <- as.data.frame(genedose)


colnames(genedose) <- genedose[1,]
genedose <- genedose[-1,]
names(genedose)[names(genedose) == 'snp'] <- 'SampleID'

genedose <- as.data.frame(sapply(genedose, as.numeric))

## split colnames
## three snps named incorrectly that will break this function, rename before reading in
t <- colnames(genedose)
t <- strsplit(t, ":")
colnames(genedose) <- sapply(t, "[", 1)

## get names of snps
snpNames<-colnames(genedose)
snpNames


## join geneDose with residuals
STAT_residuals_geneDose <- merge(STAT_residuals, genedose, by = "SampleID")
View(STAT_residuals_geneDose)



model_stats<-matrix(ncol=6,nrow=length(snpNames))

for (i in 2:length(snpNames)) {
  print(snpNames[i])
}
  
for (i in 2:length(snpNames)) {
  formula<-as.formula( paste(snpNames[i], paste( "residual" ), sep=" ~ " ) )
  print(formula)
  model<-glm(formula, data=STAT_residuals_geneDose)
  
  model_stats[i,1]<-as.numeric(coef(model)[2])#beta
  model_stats[i,2:3]<-as.numeric(confint(model,"residual"))#upper-, lower-CI
  model_stats[i,4]<-(model_stats[i,3] - model_stats[i,2])/(2*1.96)#SE
  model_stats[i,5]<-1/(model_stats[i,4]^2)#weights
  model_stats[i,6]<-summary(model)$coeff[2,4]#p-value
  
  
  
}

colnames(model_stats)<-c("Beta","Upper CI","Lower CI","SE","Weights","P-Value")
rownames(model_stats)<-snpNames	
View(model_stats)
write.csv(model_stats, "C:\\Users\\alan_\\Desktop\\RAanalysis\\lungAcuteResidualsSnp.csv")
summary(model_stats)




###### associate prs and wprs with toxicity endpoints
toxEndPoints <- colnames(toxicityFiltered)
toxEndPoints <- toxEndPoints[-c(1, 2)]
toxEndPoints




model_stats_toxPRS<-matrix(ncol=6,nrow=(length(toxEndPoints)))

for (i in 1:length(toxEndPoints)) {
  tmp <- paste("factor(",toxEndPoints[i],")")
  formula<-as.formula( paste(tmp, paste( "prs_alan" ), sep=" ~ " ) )
  print(formula)
  model<-glm(formula, family=binomial(link='logit'), data=STAT_residuals_geneDose)
  #print(model)
  model_stats_toxPRS[i,1]<-as.numeric(coef(model)[2])#beta
  #print(model_stats_toxPRS[i,1])
  model_stats_toxPRS[i,2:3]<-as.numeric(confint(model,"prs_alan"))#upper-, lower-CI
  #print(model_stats_toxPRS[i,2:3])
  model_stats_toxPRS[i,4]<-(model_stats[i,3] - model_stats[i,2])/(2*1.96)#SE
  model_stats_toxPRS[i,5]<-1/(model_stats[i,4]^2)#weights
  model_stats_toxPRS[i,6]<-summary(model)$coeff[2,4]#p-value
  
  
}

colnames(model_stats_toxPRS)<-c("Beta","Upper CI","Lower CI","SE","Weights","P-Value")
rownames(model_stats_toxPRS)<-toxEndPoints	
View(model_stats_toxPRS)
write.csv(model_stats_toxPRS, "C:\\Users\\alan_\\Desktop\\RAanalysis\\prostateAcuteResults_toxictyPRS.csv")
#summary(model_stats_toxPRS)

###################

model_stats_toxPRS<-matrix(ncol=6,nrow=(length(toxEndPoints)))

for (i in 1:length(toxEndPoints)) {
  tmp <- paste("factor(",toxEndPoints[i],")")
  formula<-as.formula( paste(tmp, paste( "prs_precentile_alan" ), sep=" ~ " ) )
  print(formula)
  model<-glm(formula, family=binomial(link='logit'), data=STAT_residuals_geneDose)
  #print(model)
  model_stats_toxPRS[i,1]<-as.numeric(coef(model)[2])#beta
  #print(model_stats_toxPRS[i,1])
  #model_stats_toxPRS[i,2:3]<-as.numeric(confint(model,"prs_precentile_alan"))#upper-, lower-CI
  #print(model_stats_toxPRS[i,2:3])
  model_stats_toxPRS[i,4]<-(model_stats[i,3] - model_stats[i,2])/(2*1.96)#SE
  model_stats_toxPRS[i,5]<-1/(model_stats[i,4]^2)#weights
  model_stats_toxPRS[i,6]<-summary(model)$coeff[2,4]#p-value
  
  
}

colnames(model_stats_toxPRS)<-c("Beta","Upper CI","Lower CI","SE","Weights","P-Value")
rownames(model_stats_toxPRS)<-toxEndPoints	
View(model_stats_toxPRS)
write.csv(model_stats_toxPRS, "C:\\Users\\alan_\\Desktop\\RAanalysis\\prostateAcuteResults_toxictyPRS_percentile.csv")



######
model_stats_toxWPRS<-matrix(ncol=6,nrow=(length(toxEndPoints)))

for (i in 1:length(toxEndPoints)) {
  tmp <- paste("factor(",toxEndPoints[i],")")
  formula<-as.formula( paste(tmp, paste( "wprs_alan" ), sep=" ~ " ) )
  print(formula)
  model<-glm(formula, family=binomial(link='logit'), data=STAT_residuals_geneDose)
  
  model_stats_toxWPRS[i,1]<-as.numeric(coef(model)[2])#beta
  #print(model_stats_toxPRS[i,1])
  model_stats_toxWPRS[i,2:3]<-as.numeric(confint(model,"wprs_alan"))#upper-, lower-CI
  #print(model_stats_toxPRS[i,2:3])
  model_stats_toxWPRS[i,4]<-(model_stats[i,3] - model_stats[i,2])/(2*1.96)#SE
  model_stats_toxWPRS[i,5]<-1/(model_stats[i,4]^2)#weights
  model_stats_toxWPRS[i,6]<-summary(model)$coeff[2,4]#p-value
  
  
}

colnames(model_stats_toxWPRS)<-c("Beta","Upper CI","Lower CI","SE","Weights","P-Value")
rownames(model_stats_toxWPRS)<-toxEndPoints	
View(model_stats_toxWPRS)
write.csv(model_stats_toxWPRS, "C:\\Users\\alan_\\Desktop\\RAanalysis\\prostateAcuteResults_toxictyWPRS.csv")
summary(model_stats_toxWPRS)


######
model_stats_toxWPRS<-matrix(ncol=6,nrow=(length(toxEndPoints)))

for (i in 1:length(toxEndPoints)) {
  tmp <- paste("factor(",toxEndPoints[i],")")
  formula<-as.formula( paste(tmp, paste( "wprs_precentile_alan" ), sep=" ~ " ) )
  print(formula)
  model<-glm(formula, family=binomial(link='logit'), data=STAT_residuals_geneDose)
  
  model_stats_toxWPRS[i,1]<-as.numeric(coef(model)[2])#beta
  #print(model_stats_toxPRS[i,1])
  #model_stats_toxWPRS[i,2:3]<-as.numeric(confint(model,"wprs_precentile_alan"))#upper-, lower-CI
  #print(model_stats_toxPRS[i,2:3])
  model_stats_toxWPRS[i,4]<-(model_stats[i,3] - model_stats[i,2])/(2*1.96)#SE
  model_stats_toxWPRS[i,5]<-1/(model_stats[i,4]^2)#weights
  model_stats_toxWPRS[i,6]<-summary(model)$coeff[2,4]#p-value
  
  
}

colnames(model_stats_toxWPRS)<-c("Beta","Upper CI","Lower CI","SE","Weights","P-Value")
rownames(model_stats_toxWPRS)<-toxEndPoints	
View(model_stats_toxWPRS)
write.csv(model_stats_toxWPRS, "C:\\Users\\alan_\\Desktop\\RAanalysis\\prostateAcuteResults_toxictyWPRS_percentile.csv")
summary(model_stats_toxWPRS)







######################################
######################################
######################################
######################################


plot(STAT_prs_factors$age_at_radiotherapy_start_yrs, STAT_prs_factors$prs)

### split of age for wprs
summary(STAT_prs_factors$age_at_radiotherapy_start_yrs)
STAT_prs_factors$ageSplit <- STAT_prs_factors$age_at_radiotherapy_start_yrs > median((STAT_prs_factors$age_at_radiotherapy_start_yrs))

tapply(STAT_prs_factors$wprs, STAT_prs_factors$ageSplit, summary)
tapply(STAT_prs_factors$prs, STAT_prs_factors$ageSplit, summary)

ggplot(STAT_prs_factors, aes(x = ageSplit, y = wprs)) + 
  geom_boxplot() + 
  theme_classic() +
  stat_compare_means(method = "wilcox.test", aes(group = ageSplit, label = paste0("p = ",..p.format..)), label.x = 1.4, label.y = 6, size = 6)

ggplot(STAT_prs_factors, aes(x = ageSplit, y = prs)) + 
  geom_boxplot() + 
  theme_classic() +
  stat_compare_means(method = "wilcox.test", aes(group = ageSplit, label = paste0("p = ",..p.format..)), label.x = 1.4, label.y = 6, size = 6)

plot(STAT_prs_factors$age_at_radiotherapy_start_yrs, STAT_prs_factors$prs)


summary(STAT_prs_factors$Country)
tapply(STAT_prs_factors$wprs, STAT_prs_factors$Country, summary)
tapply(STAT_prs_factors$STAT, STAT_prs_factors$Country, summary)


ggplot(STAT_prs_factors, aes(x = Country, y = wprs)) + 
  geom_boxplot() + 
  theme_classic() +
  stat_compare_means(method = "wilcox.test", aes(group = Country, label = paste0("p = ",..p.format..)), label.x = 1.4, label.y = 6, size = 6)

ggplot(STAT_prs_factors, aes(x = Country, y = prs)) + 
  geom_boxplot() + 
  theme_classic() +
  stat_compare_means(method = "wilcox.test", aes(group = Country, label = paste0("p = ",..p.format..)), label.x = 1.4, label.y = 6, size = 6)




#############################
sarahList <- read.csv("C:\\Users\\alan_\\Desktop\\RAanalysis\\listSarah\\REQUITE_prostate_STATacute_for_Alan.txt", sep = '\t')
View(sarahList)
colnames(sarahList)
colnames(sarahList) <- c("SubjectId", "array_id", "stat_acute", "rstat_acute")
colnames(sarahList)
colnames(STAT_prs_factors)
View(STAT_prs_factors)
STATcomp <- merge(STAT_prs_factors, sarahList, by = 'SubjectId')
STATcomp <- merge(toxicityFilteredSTAT, sarahList, by = 'SubjectId')

View(STATcomp)

STATcomp$STATdiff <- STATcomp$STAT - STATcomp$stat_acute
summary(STATcomp$STATdiff)
View(STATcomp)

ttt <- glm(STAT~wprs, data = STATcomp)
summary(ttt)
ttt <- glm(stat_acute~wprs, data = STATcomp)
summary(ttt)

ttt <- glm(STAT~wprs + age_at_radiotherapy_start_yrs + diabetes + p3radical_prostatectomy + p3hormone_therapy + doseBED + ra, data = STATcomp)
summary(ttt) 
ttt <- glm(stat_acute~wprs + age_at_radiotherapy_start_yrs + diabetes + p3radical_prostatectomy + p3hormone_therapy + doseBED + ra, data = STATcomp)
summary(ttt) 
