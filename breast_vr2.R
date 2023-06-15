

############################################################################################
####
#### A McWilliam
#### vr2 27 Feb 2023 updated for new breast analysis
#### 
#### PRO
####
############################################################################################


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
Tox <- read.csv("C:/Users/alan_/Desktop/RAanalysis/REQUITEdata/Breast/study/datasets/dataset5016.tsv", sep="\t", header=T )
View(Tox)
Treat <- read.csv("C:/Users/alan_/Desktop/RAanalysis/REQUITEdata/Breast/study/datasets/dataset5018.tsv", sep="\t", header=T)
View(Treat)
Factor <- read.csv("C:/Users/alan_/Desktop/RAanalysis/REQUITEdata/Breast/study/datasets/dataset5021.tsv", sep="\t", header=T)
View(Factor)

#lungPRO <- read.csv("C:/Users/alan_/Desktop/RAanalysis/REQUITEdata/Lung/study/datasets/dataset5004.tsv", sep="\t", header=T)
#View(lungPRO)

### select ID and radiotherapy start data
radiotherapyStart <- Treat %>%
  select(SubjectId, b3radio_breast_startdate)


### select out toxicity data of interest for STAT
## breast_erythema, breast_skin_ulceration -- ACUTE
#Telangiectasia, Oedema, atrophy, Induration,shrinkage Pigmentation, 
toxicity <- Tox %>%
  select(SubjectId, date, breast_telangiectasia_outside_tumour_bed, breast_telangiectasia_tumour_bed, breast_skin_induration_outside_tumour_bed, breast_skin_induration_tumour_bed, breast_skin_hyperpigmentation, breast_atrophy, breast_oedema)

##breast_telangiectasia_tumour_bed, breast_skin_induration_tumour_bed

toxicity <- merge(Tox, radiotherapyStart, by = 'SubjectId')
toxicity$monthStartTreat <- interval(ymd(as.Date(toxicity$b3radio_breast_startdate)), ymd(as.Date(toxicity$date)))
toxicity$monthStartTreat = toxicity$monthStartTreat %/% months(1)
### 

### need to set a time point for selecting highest toxicity score and filter
### keep highest recorded score for each toxicity 
months = 500 ## acute toxicity - 90 days acute set 
minMonths = 3

## identify patients with baseline values
baselineCounts <- toxicity %>%
  group_by(SubjectId) %>%
  filter(monthStartTreat == 0) %>%
  select(SubjectId)


## filter by baseline and subtract first entry
toxicitySubtract <- toxicity %>%
  select(SubjectId, monthStartTreat, breast_telangiectasia_outside_tumour_bed, breast_telangiectasia_tumour_bed, breast_skin_induration_outside_tumour_bed, breast_skin_induration_tumour_bed, breast_skin_hyperpigmentation, breast_atrophy, breast_oedema) %>%
  filter(SubjectId %in% baselineCounts$SubjectId) %>%
  group_by(SubjectId) %>%
  mutate(telangiectasia_diff = breast_telangiectasia_outside_tumour_bed - first(breast_telangiectasia_outside_tumour_bed), 
         telangiectasiaTumour_diff = breast_telangiectasia_tumour_bed - first(breast_telangiectasia_tumour_bed),
         induration_diff = breast_skin_induration_outside_tumour_bed - first(breast_skin_induration_outside_tumour_bed),
         indurationTumour_diff = breast_skin_induration_tumour_bed - first(breast_skin_induration_tumour_bed),
         pigmentation_diff = breast_skin_hyperpigmentation - first(breast_skin_hyperpigmentation), 
         atrophy_diff = breast_atrophy - first(breast_atrophy),
         oedema_diff = breast_oedema - first(breast_oedema)) %>%
  filter(monthStartTreat > minMonths) %>%
  filter(monthStartTreat < months) %>%
  summarise(maxMonths = max(monthStartTreat), 
            maxTelangiectasia = max(telangiectasia_diff),
            maxTelangiectasiaTumour = max(telangiectasiaTumour_diff),
            maxInduration = max(induration_diff), 
            maxIndurationTumour = max(indurationTumour_diff),
            maxPigmentation = max(pigmentation_diff), 
            maxAtrophy = max(atrophy_diff),
            maxOedema = max(oedema_diff)) %>%
  select(SubjectId, maxMonths, maxTelangiectasia, maxTelangiectasiaTumour, maxInduration, maxIndurationTumour, maxPigmentation, maxAtrophy, maxOedema)

toxicitySubtract[toxicitySubtract < 0] <- 0 
View(toxicitySubtract)


## patients without baseline
toxicityNoBaseline <- toxicity %>%
  filter(!SubjectId %in% baselineCounts$SubjectId) %>%
  filter(monthStartTreat > minMonths) %>%
  filter(monthStartTreat < months) %>%
  group_by(SubjectId) %>%
  summarise(maxMonths = max(monthStartTreat), 
            maxTelangiectasia = max(breast_telangiectasia_outside_tumour_bed),
            maxTelangiectasiaTumour = max(breast_telangiectasia_tumour_bed),
            maxInduration = max(breast_skin_induration_outside_tumour_bed),
            maxIndurationTumour = max(breast_skin_induration_tumour_bed),
            maxPigmentation = max(breast_skin_hyperpigmentation), 
            maxAtrophy = max(breast_atrophy),
            maxOedema = max(breast_oedema)) %>%
  select(SubjectId, maxMonths, maxTelangiectasia, maxTelangiectasiaTumour, maxInduration, maxIndurationTumour, maxPigmentation, maxAtrophy, maxOedema)


## join both together again
toxicityFiltered <- rbind(toxicitySubtract, toxicityNoBaseline)
View(toxicityFiltered)


###create temp data frame to display summary stats
toxicity_summaryStats <- toxicityFiltered %>%
  select(maxTelangiectasia, maxTelangiectasiaTumour, maxInduration, maxIndurationTumour, maxPigmentation, maxAtrophy, maxOedema)

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

#toxicityFilteredSTAT <- toxicityFilteredSTAT[!is.nan(toxicityFilteredSTAT$STAT), ]
View(toxicityFilteredSTAT)
write.csv(toxicityFiltered, "C:/Users/alan_/Desktop/breastLate.csv")

### select other clinical variables for inclusion in analysis


### save csv of data wkith STAT ready for analysis.
#write.csv(toxicityFilteredSTAT, "C:/Users/alan_/Desktop/RAanalysis/REQUITEdata/processed/prostate_acuteSTAT.csv")


############################################################################################
### summary of patient characteristics
### summary STAT and plot histogram
###

summary(toxicityFilteredSTAT$STAT)

ggplot(data=toxicityFilteredSTAT, aes(STAT)) + 
  geom_histogram(breaks=seq(-1,3, by = 0.25),
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
sampleIDlink <- read.csv("C:/Users/alan_/Desktop/RAanalysis/REQUITEdata/Breast/study/datasets/dataset5039.tsv", sep="\t", header=T)

##t <- merge(sampleIDlink, toxicityCROFilteredSTAT, by = "SubjectId")
STAT_prs <- merge(PRS_all, sampleIDlink, by = "SampleID")
STAT_prs <- merge(STAT_prs, toxicityFilteredSTAT, by = "SubjectId")
View(STAT_prs)


#t <- glm(STAT~prs_precentile_alan, data = STAT_prs)
#STAT_prs$prs_precentile_alan <- STAT_prs$prs_alan > quantile(STAT_prs$prs_alan, c(.95)) 
#summary(t)

##########################################################################################
### need to select other patient factors
### BED for acute tox - 3 also try 10 for normal tissue???

##Age
##Smoker - Y/N
##Chemotherapy - Y/N
##Cardiovascular disease - Y/N
##BMI
##Breast volume (cm3)
##Diabetes - Y/N
##Post op breast infection -Y/N
##Breast boost - Y/N

patInfo <- read.csv("C:/Users/alan_/Desktop/RAanalysis/REQUITEdata/Breast/study/datasets/dataset5020.tsv", sep="\t", header=T)
#View(patInfo)

patTreat <- Treat %>%
  select( SubjectId, b3radio_breast_fractions, b3radio_breast_dose_Gy, b3chemo_adjuvant, b3chemo_neo_adjuvant, b3radio_boost, b3post_operative_infection, b3radio_breast_ct_volume_cm3 ) 

patDetail <- patInfo %>%
  select(SubjectId, age_at_radiotherapy_start_yrs, smoker, height_cm, weight_at_cancer_diagnosis_kg, history_of_heart_disease, ra, diabetes)

## calc BMI
patDetail$BMI <- patDetail$weight_at_cancer_diagnosis_kg / ((patDetail$height_cm / 100))^2 


### need to calculate the BED prescribed
### BED = D x (1 + [d / (??/??)])
alphaBeta <- 10 #3
patTreat$doseFraction <- patTreat$b3radio_breast_dose_Gy/patTreat$b3radio_breast_fractions
patTreat$doseBED <-   patTreat$b3radio_breast_dose_Gy  * (1 + (patTreat$doseFraction/alphaBeta))

patDetail$smoker[patDetail$smoker > 0 ] <- 1 

#STAT_prs_factors <- merge(STAT_prs, patFactors, by = "SubjectId")
STAT_prs_factors <- merge(STAT_prs, patTreat, by = "SubjectId")
STAT_prs_factors <- merge(STAT_prs_factors, patDetail, by = "SubjectId")

STAT_prs_factors$chemo = STAT_prs_factors$b3chemo_adjuvant + STAT_prs_factors$b3chemo_neo_adjuvant
STAT_prs_factors$chemo[STAT_prs_factors$chemo > 0 ] <- 1



View(STAT_prs_factors)

view(dfSummary(STAT_prs_factors))



##########################################################################################
### analysing PRS and wPRS

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




t <- glm(STAT~prs_alan + age_at_radiotherapy_start_yrs + factor(smoker) + factor(chemo) + factor(history_of_heart_disease) + BMI + b3radio_breast_ct_volume_cm3 + doseBED +  factor(b3radio_boost) + factor(b3post_operative_infection) + factor(diabetes) + factor(ra), data = STAT_prs_factors)
summary(t) 
confint(t)
t <- glm(STAT~wprs_alan + age_at_radiotherapy_start_yrs + factor(smoker) + factor(chemo) + factor(history_of_heart_disease) + BMI + b3radio_breast_ct_volume_cm3 + doseBED +  factor(b3radio_boost) + factor(b3post_operative_infection) + factor(diabetes) + factor(ra), data = STAT_prs_factors)
summary(t) 
confint(t)

t <- glm(STAT~prs_precentile_alan  + age_at_radiotherapy_start_yrs + factor(smoker) + factor(chemo) + factor(history_of_heart_disease) + BMI + b3radio_breast_ct_volume_cm3 + doseBED +  factor(b3radio_boost) + factor(b3post_operative_infection) + factor(diabetes) + factor(ra), data = STAT_prs_factors)
summary(t) 
confint(t)
t <- glm(STAT~wprs_precentile_alan + age_at_radiotherapy_start_yrs + factor(smoker) + factor(chemo) + factor(history_of_heart_disease) + BMI + b3radio_breast_ct_volume_cm3 + doseBED +  factor(b3radio_boost) + factor(b3post_operative_infection) + factor(diabetes) + factor(ra), data = STAT_prs_factors)
summary(t) 
confint(t)


##### individual endpoints

## maxTelangiectasia, 
#prs
tt1 <- glm(maxTelangiectasia~prs_alan + age_at_radiotherapy_start_yrs + factor(smoker) + factor(chemo) + factor(history_of_heart_disease) + BMI + b3radio_breast_ct_volume_cm3 + doseBED +  factor(b3radio_boost) + factor(b3post_operative_infection) + factor(diabetes) + factor(ra), data = STAT_prs_factors)
summary(tt1) 
ttt1 <- confint(tt1, level = 0.95)
a1 <- paste(format(round(ttt1[2,1], 3), nsmall = 3), format(round(ttt1[2,2], 3), nsmall = 3), sep=', ')
aa1<- paste('(',a1,')', sep = '')
aaa1 <- paste(format(round(as.numeric(summary(tt1)$coeff[2]), 3), nsmall = 3), aa1, format(round(summary(tt1)$coeff[2,4], 2), nsmall = 2), sep = ' ')


#wprs
tt1 <-glm(maxTelangiectasia~wprs_alan + age_at_radiotherapy_start_yrs + factor(smoker) + factor(chemo) + factor(history_of_heart_disease) + BMI + b3radio_breast_ct_volume_cm3 + doseBED +  factor(b3radio_boost) + factor(b3post_operative_infection) + factor(diabetes) + factor(ra), data = STAT_prs_factors)
summary(tt1) 
ttt1 <- confint(tt1, level = 0.95)
a2 <- paste(format(round(ttt1[2,1], 3), nsmall = 3), format(round(ttt1[2,2], 3), nsmall = 3), sep=', ')
aa2<- paste('(',a2,')', sep = '')
aaa2 <- paste(format(round(as.numeric(summary(tt1)$coeff[2]), 3), nsmall = 3), aa2, format(round(summary(tt1)$coeff[2,4], 2), nsmall = 2), sep = ' ')

aaa1
aaa2

#prs percentile
tt1 <- glm(maxTelangiectasia~prs_precentile_alan + age_at_radiotherapy_start_yrs + factor(smoker) + factor(chemo) + factor(history_of_heart_disease) + BMI + b3radio_breast_ct_volume_cm3 + doseBED +  factor(b3radio_boost) + factor(b3post_operative_infection) + factor(diabetes) + factor(ra), data = STAT_prs_factors)
summary(tt1) 
ttt1 <- confint(tt1, level = 0.95)
a1 <- paste(format(round(ttt1[2,1], 3), nsmall = 3), format(round(ttt1[2,2], 3), nsmall = 3), sep=', ')
aa1<- paste('(',a1,')', sep = '')
aaa1 <- paste(format(round(as.numeric(summary(tt1)$coeff[2]), 3), nsmall = 3), aa1, format(round(summary(tt1)$coeff[2,4], 2), nsmall = 2), sep = ' ')


#wprs percentile
tt1 <- glm(maxTelangiectasia~wprs_precentile_alan + age_at_radiotherapy_start_yrs + factor(smoker) + factor(chemo) + factor(history_of_heart_disease) + BMI + b3radio_breast_ct_volume_cm3 + doseBED +  factor(b3radio_boost) + factor(b3post_operative_infection) + factor(diabetes) + factor(ra), data = STAT_prs_factors)
summary(tt1) 
ttt1 <- confint(tt1, level = 0.95)
a2 <- paste(format(round(ttt1[2,1], 3), nsmall = 3), format(round(ttt1[2,2], 3), nsmall = 3), sep=', ')
aa2<- paste('(',a2,')', sep = '')
aaa2 <- paste(format(round(as.numeric(summary(tt1)$coeff[2]), 3), nsmall = 3), aa2, format(round(summary(tt1)$coeff[2,4], 2), nsmall = 2), sep = ' ')

aaa1
aaa2

## maxInduration, 
#prs
tt1 <- glm(maxInduration~prs_alan + age_at_radiotherapy_start_yrs + factor(smoker) + factor(chemo) + factor(history_of_heart_disease) + BMI + b3radio_breast_ct_volume_cm3 + doseBED +  factor(b3radio_boost) + factor(b3post_operative_infection) + factor(diabetes) + factor(ra), data = STAT_prs_factors)
summary(tt1) 
ttt1 <- confint(tt1, level = 0.95)
a1 <- paste(format(round(ttt1[2,1], 3), nsmall = 3), format(round(ttt1[2,2], 3), nsmall = 3), sep=', ')
aa1<- paste('(',a1,')', sep = '')
aaa1 <- paste(format(round(as.numeric(summary(tt1)$coeff[2]), 3), nsmall = 3), aa1, format(round(summary(tt1)$coeff[2,4], 2), nsmall = 2), sep = ' ')


#wprs
tt1 <-glm(maxInduration~wprs_alan + age_at_radiotherapy_start_yrs + factor(smoker) + factor(chemo) + factor(history_of_heart_disease) + BMI + b3radio_breast_ct_volume_cm3 + doseBED +  factor(b3radio_boost) + factor(b3post_operative_infection) + factor(diabetes) + factor(ra), data = STAT_prs_factors)
summary(tt1) 
ttt1 <- confint(tt1, level = 0.95)
a2 <- paste(format(round(ttt1[2,1], 3), nsmall = 3), format(round(ttt1[2,2], 3), nsmall = 3), sep=', ')
aa2<- paste('(',a2,')', sep = '')
aaa2 <- paste(format(round(as.numeric(summary(tt1)$coeff[2]), 3), nsmall = 3), aa2, format(round(summary(tt1)$coeff[2,4], 2), nsmall = 2), sep = ' ')

aaa1
aaa2

#prs percentile
tt1 <- glm(maxInduration~prs_precentile_alan + age_at_radiotherapy_start_yrs + factor(smoker) + factor(chemo) + factor(history_of_heart_disease) + BMI + b3radio_breast_ct_volume_cm3 + doseBED +  factor(b3radio_boost) + factor(b3post_operative_infection) + factor(diabetes) + factor(ra), data = STAT_prs_factors)
summary(tt1) 
ttt1 <- confint(tt1, level = 0.95)
a1 <- paste(format(round(ttt1[2,1], 3), nsmall = 3), format(round(ttt1[2,2], 3), nsmall = 3), sep=', ')
aa1<- paste('(',a1,')', sep = '')
aaa1 <- paste(format(round(as.numeric(summary(tt1)$coeff[2]), 3), nsmall = 3), aa1, format(round(summary(tt1)$coeff[2,4], 2), nsmall = 2), sep = ' ')


#wprs percentile
tt1 <- glm(maxInduration~wprs_precentile_alan + age_at_radiotherapy_start_yrs + factor(smoker) + factor(chemo) + factor(history_of_heart_disease) + BMI + b3radio_breast_ct_volume_cm3 + doseBED +  factor(b3radio_boost) + factor(b3post_operative_infection) + factor(diabetes) + factor(ra), data = STAT_prs_factors)
summary(tt1) 
ttt1 <- confint(tt1, level = 0.95)
a2 <- paste(format(round(ttt1[2,1], 3), nsmall = 3), format(round(ttt1[2,2], 3), nsmall = 3), sep=', ')
aa2<- paste('(',a2,')', sep = '')
aaa2 <- paste(format(round(as.numeric(summary(tt1)$coeff[2]), 3), nsmall = 3), aa2, format(round(summary(tt1)$coeff[2,4], 2), nsmall = 2), sep = ' ')

aaa1
aaa2

##maxPigmentation, 
#prs
tt1 <- glm(maxPigmentation~prs_alan + age_at_radiotherapy_start_yrs + factor(smoker) + factor(chemo) + factor(history_of_heart_disease) + BMI + b3radio_breast_ct_volume_cm3 + doseBED +  factor(b3radio_boost) + factor(b3post_operative_infection) + factor(diabetes) + factor(ra), data = STAT_prs_factors)
summary(tt1) 
ttt1 <- confint(tt1, level = 0.95)
a1 <- paste(format(round(ttt1[2,1], 3), nsmall = 3), format(round(ttt1[2,2], 3), nsmall = 3), sep=', ')
aa1<- paste('(',a1,')', sep = '')
aaa1 <- paste(format(round(as.numeric(summary(tt1)$coeff[2]), 3), nsmall = 3), aa1, format(round(summary(tt1)$coeff[2,4], 2), nsmall = 2), sep = ' ')


#wprs
tt1 <-glm(maxPigmentation~wprs_alan + age_at_radiotherapy_start_yrs + factor(smoker) + factor(chemo) + factor(history_of_heart_disease) + BMI + b3radio_breast_ct_volume_cm3 + doseBED +  factor(b3radio_boost) + factor(b3post_operative_infection) + factor(diabetes) + factor(ra), data = STAT_prs_factors)
summary(tt1) 
ttt1 <- confint(tt1, level = 0.95)
a2 <- paste(format(round(ttt1[2,1], 3), nsmall = 3), format(round(ttt1[2,2], 3), nsmall = 3), sep=', ')
aa2<- paste('(',a2,')', sep = '')
aaa2 <- paste(format(round(as.numeric(summary(tt1)$coeff[2]), 3), nsmall = 3), aa2, format(round(summary(tt1)$coeff[2,4], 2), nsmall = 2), sep = ' ')

aaa1
aaa2

#prs percentile
tt1 <- glm(maxPigmentation~prs_precentile_alan + age_at_radiotherapy_start_yrs + factor(smoker) + factor(chemo) + factor(history_of_heart_disease) + BMI + b3radio_breast_ct_volume_cm3 + doseBED +  factor(b3radio_boost) + factor(b3post_operative_infection) + factor(diabetes) + factor(ra), data = STAT_prs_factors)
summary(tt1) 
ttt1 <- confint(tt1, level = 0.95)
a1 <- paste(format(round(ttt1[2,1], 3), nsmall = 3), format(round(ttt1[2,2], 3), nsmall = 3), sep=', ')
aa1<- paste('(',a1,')', sep = '')
aaa1 <- paste(format(round(as.numeric(summary(tt1)$coeff[2]), 3), nsmall = 3), aa1, format(round(summary(tt1)$coeff[2,4], 2), nsmall = 2), sep = ' ')


#wprs percentile
tt1 <- glm(maxPigmentation~wprs_precentile_alan + age_at_radiotherapy_start_yrs + factor(smoker) + factor(chemo) + factor(history_of_heart_disease) + BMI + b3radio_breast_ct_volume_cm3 + doseBED +  factor(b3radio_boost) + factor(b3post_operative_infection) + factor(diabetes) + factor(ra), data = STAT_prs_factors)
summary(tt1) 
ttt1 <- confint(tt1, level = 0.95)
a2 <- paste(format(round(ttt1[2,1], 3), nsmall = 3), format(round(ttt1[2,2], 3), nsmall = 3), sep=', ')
aa2<- paste('(',a2,')', sep = '')
aaa2 <- paste(format(round(as.numeric(summary(tt1)$coeff[2]), 3), nsmall = 3), aa2, format(round(summary(tt1)$coeff[2,4], 2), nsmall = 2), sep = ' ')

aaa1
aaa2

##maxAtrophy, 
#prs
tt1 <- glm(maxAtrophy~prs_alan + age_at_radiotherapy_start_yrs + factor(smoker) + factor(chemo) + factor(history_of_heart_disease) + BMI + b3radio_breast_ct_volume_cm3 + doseBED +  factor(b3radio_boost) + factor(b3post_operative_infection) + factor(diabetes) + factor(ra), data = STAT_prs_factors)
summary(tt1) 
ttt1 <- confint(tt1, level = 0.95)
a1 <- paste(format(round(ttt1[2,1], 3), nsmall = 3), format(round(ttt1[2,2], 3), nsmall = 3), sep=', ')
aa1<- paste('(',a1,')', sep = '')
aaa1 <- paste(format(round(as.numeric(summary(tt1)$coeff[2]), 3), nsmall = 3), aa1, format(round(summary(tt1)$coeff[2,4], 2), nsmall = 2), sep = ' ')


#wprs
tt1 <-glm(maxAtrophy~wprs_alan + age_at_radiotherapy_start_yrs + factor(smoker) + factor(chemo) + factor(history_of_heart_disease) + BMI + b3radio_breast_ct_volume_cm3 + doseBED +  factor(b3radio_boost) + factor(b3post_operative_infection) + factor(diabetes) + factor(ra), data = STAT_prs_factors)
summary(tt1) 
ttt1 <- confint(tt1, level = 0.95)
a2 <- paste(format(round(ttt1[2,1], 3), nsmall = 3), format(round(ttt1[2,2], 3), nsmall = 3), sep=', ')
aa2<- paste('(',a2,')', sep = '')
aaa2 <- paste(format(round(as.numeric(summary(tt1)$coeff[2]), 3), nsmall = 3), aa2, format(round(summary(tt1)$coeff[2,4], 2), nsmall = 2), sep = ' ')

aaa1
aaa2

#prs percentile
tt1 <- glm(maxAtrophy~prs_precentile_alan + age_at_radiotherapy_start_yrs + factor(smoker) + factor(chemo) + factor(history_of_heart_disease) + BMI + b3radio_breast_ct_volume_cm3 + doseBED +  factor(b3radio_boost) + factor(b3post_operative_infection) + factor(diabetes) + factor(ra), data = STAT_prs_factors)
summary(tt1) 
ttt1 <- confint(tt1, level = 0.95)
a1 <- paste(format(round(ttt1[2,1], 3), nsmall = 3), format(round(ttt1[2,2], 3), nsmall = 3), sep=', ')
aa1<- paste('(',a1,')', sep = '')
aaa1 <- paste(format(round(as.numeric(summary(tt1)$coeff[2]), 3), nsmall = 3), aa1, format(round(summary(tt1)$coeff[2,4], 2), nsmall = 2), sep = ' ')


#wprs percentile
tt1 <- glm(maxAtrophy~wprs_precentile_alan + age_at_radiotherapy_start_yrs + factor(smoker) + factor(chemo) + factor(history_of_heart_disease) + BMI + b3radio_breast_ct_volume_cm3 + doseBED +  factor(b3radio_boost) + factor(b3post_operative_infection) + factor(diabetes) + factor(ra), data = STAT_prs_factors)
summary(tt1) 
ttt1 <- confint(tt1, level = 0.95)
a2 <- paste(format(round(ttt1[2,1], 3), nsmall = 3), format(round(ttt1[2,2], 3), nsmall = 3), sep=', ')
aa2<- paste('(',a2,')', sep = '')
aaa2 <- paste(format(round(as.numeric(summary(tt1)$coeff[2]), 3), nsmall = 3), aa2, format(round(summary(tt1)$coeff[2,4], 2), nsmall = 2), sep = ' ')

aaa1
aaa2

## maxOedema

#prs
tt1 <- glm(maxOedema~prs_alan + age_at_radiotherapy_start_yrs + factor(smoker) + factor(chemo) + factor(history_of_heart_disease) + BMI + b3radio_breast_ct_volume_cm3 + doseBED +  factor(b3radio_boost) + factor(b3post_operative_infection) + factor(diabetes) + factor(ra), data = STAT_prs_factors)
summary(tt1) 
ttt1 <- confint(tt1, level = 0.95)
a1 <- paste(format(round(ttt1[2,1], 3), nsmall = 3), format(round(ttt1[2,2], 3), nsmall = 3), sep=', ')
aa1<- paste('(',a1,')', sep = '')
aaa1 <- paste(format(round(as.numeric(summary(tt1)$coeff[2]), 3), nsmall = 3), aa1, format(round(summary(tt1)$coeff[2,4], 2), nsmall = 2), sep = ' ')


#wprs
tt1 <-glm(maxOedema~wprs_alan + age_at_radiotherapy_start_yrs + factor(smoker) + factor(chemo) + factor(history_of_heart_disease) + BMI + b3radio_breast_ct_volume_cm3 + doseBED +  factor(b3radio_boost) + factor(b3post_operative_infection) + factor(diabetes) + factor(ra), data = STAT_prs_factors)
summary(tt1) 
ttt1 <- confint(tt1, level = 0.95)
a2 <- paste(format(round(ttt1[2,1], 3), nsmall = 3), format(round(ttt1[2,2], 3), nsmall = 3), sep=', ')
aa2<- paste('(',a2,')', sep = '')
aaa2 <- paste(format(round(as.numeric(summary(tt1)$coeff[2]), 3), nsmall = 3), aa2, format(round(summary(tt1)$coeff[2,4], 2), nsmall = 2), sep = ' ')

aaa1
aaa2

#prs percentile
tt1 <- glm(maxOedema~prs_precentile_alan + age_at_radiotherapy_start_yrs + factor(smoker) + factor(chemo) + factor(history_of_heart_disease) + BMI + b3radio_breast_ct_volume_cm3 + doseBED +  factor(b3radio_boost) + factor(b3post_operative_infection) + factor(diabetes) + factor(ra), data = STAT_prs_factors)
summary(tt1) 
ttt1 <- confint(tt1, level = 0.95)
a1 <- paste(format(round(ttt1[2,1], 3), nsmall = 3), format(round(ttt1[2,2], 3), nsmall = 3), sep=', ')
aa1<- paste('(',a1,')', sep = '')
aaa1 <- paste(format(round(as.numeric(summary(tt1)$coeff[2]), 3), nsmall = 3), aa1, format(round(summary(tt1)$coeff[2,4], 2), nsmall = 2), sep = ' ')


#wprs percentile
tt1 <- glm(maxOedema~wprs_precentile_alan + age_at_radiotherapy_start_yrs + factor(smoker) + factor(chemo) + factor(history_of_heart_disease) + BMI + b3radio_breast_ct_volume_cm3 + doseBED +  factor(b3radio_boost) + factor(b3post_operative_infection) + factor(diabetes) + factor(ra), data = STAT_prs_factors)
summary(tt1) 
ttt1 <- confint(tt1, level = 0.95)
a2 <- paste(format(round(ttt1[2,1], 3), nsmall = 3), format(round(ttt1[2,2], 3), nsmall = 3), sep=', ')
aa2<- paste('(',a2,')', sep = '')
aaa2 <- paste(format(round(as.numeric(summary(tt1)$coeff[2]), 3), nsmall = 3), aa2, format(round(summary(tt1)$coeff[2,4], 2), nsmall = 2), sep = ' ')

aaa1
aaa2




## maxTelangiectasiaTumour

#prs
tt1 <- glm(maxTelangiectasiaTumour~prs_alan + age_at_radiotherapy_start_yrs + factor(smoker) + factor(chemo) + factor(history_of_heart_disease) + BMI + b3radio_breast_ct_volume_cm3 + doseBED +  factor(b3radio_boost) + factor(b3post_operative_infection) + factor(diabetes) + factor(ra), data = STAT_prs_factors)
summary(tt1) 
ttt1 <- confint(tt1, level = 0.95)
a1 <- paste(format(round(ttt1[2,1], 3), nsmall = 3), format(round(ttt1[2,2], 3), nsmall = 3), sep=', ')
aa1<- paste('(',a1,')', sep = '')
aaa1 <- paste(format(round(as.numeric(summary(tt1)$coeff[2]), 3), nsmall = 3), aa1, format(round(summary(tt1)$coeff[2,4], 2), nsmall = 2), sep = ' ')


#wprs
tt1 <-glm(maxTelangiectasiaTumour~wprs_alan + age_at_radiotherapy_start_yrs + factor(smoker) + factor(chemo) + factor(history_of_heart_disease) + BMI + b3radio_breast_ct_volume_cm3 + doseBED +  factor(b3radio_boost) + factor(b3post_operative_infection) + factor(diabetes) + factor(ra), data = STAT_prs_factors)
summary(tt1) 
ttt1 <- confint(tt1, level = 0.95)
a2 <- paste(format(round(ttt1[2,1], 3), nsmall = 3), format(round(ttt1[2,2], 3), nsmall = 3), sep=', ')
aa2<- paste('(',a2,')', sep = '')
aaa2 <- paste(format(round(as.numeric(summary(tt1)$coeff[2]), 3), nsmall = 3), aa2, format(round(summary(tt1)$coeff[2,4], 2), nsmall = 2), sep = ' ')

aaa1
aaa2

#prs percentile
tt1 <- glm(maxTelangiectasiaTumour~prs_precentile_alan + age_at_radiotherapy_start_yrs + factor(smoker) + factor(chemo) + factor(history_of_heart_disease) + BMI + b3radio_breast_ct_volume_cm3 + doseBED +  factor(b3radio_boost) + factor(b3post_operative_infection) + factor(diabetes) + factor(ra), data = STAT_prs_factors)
summary(tt1) 
ttt1 <- confint(tt1, level = 0.95)
a1 <- paste(format(round(ttt1[2,1], 3), nsmall = 3), format(round(ttt1[2,2], 3), nsmall = 3), sep=', ')
aa1<- paste('(',a1,')', sep = '')
aaa1 <- paste(format(round(as.numeric(summary(tt1)$coeff[2]), 3), nsmall = 3), aa1, format(round(summary(tt1)$coeff[2,4], 2), nsmall = 2), sep = ' ')


#wprs percentile
tt1 <- glm(maxTelangiectasiaTumour~wprs_precentile_alan + age_at_radiotherapy_start_yrs + factor(smoker) + factor(chemo) + factor(history_of_heart_disease) + BMI + b3radio_breast_ct_volume_cm3 + doseBED +  factor(b3radio_boost) + factor(b3post_operative_infection) + factor(diabetes) + factor(ra), data = STAT_prs_factors)
summary(tt1) 
ttt1 <- confint(tt1, level = 0.95)
a2 <- paste(format(round(ttt1[2,1], 3), nsmall = 3), format(round(ttt1[2,2], 3), nsmall = 3), sep=', ')
aa2<- paste('(',a2,')', sep = '')
aaa2 <- paste(format(round(as.numeric(summary(tt1)$coeff[2]), 3), nsmall = 3), aa2, format(round(summary(tt1)$coeff[2,4], 2), nsmall = 2), sep = ' ')

aaa1
aaa2



## maxIndurationTumour

#prs
tt1 <- glm(maxIndurationTumour~prs_alan + age_at_radiotherapy_start_yrs + factor(smoker) + factor(chemo) + factor(history_of_heart_disease) + BMI + b3radio_breast_ct_volume_cm3 + doseBED +  factor(b3radio_boost) + factor(b3post_operative_infection) + factor(diabetes) + factor(ra), data = STAT_prs_factors)
summary(tt1) 
ttt1 <- confint(tt1, level = 0.95)
a1 <- paste(format(round(ttt1[2,1], 3), nsmall = 3), format(round(ttt1[2,2], 3), nsmall = 3), sep=', ')
aa1<- paste('(',a1,')', sep = '')
aaa1 <- paste(format(round(as.numeric(summary(tt1)$coeff[2]), 3), nsmall = 3), aa1, format(round(summary(tt1)$coeff[2,4], 2), nsmall = 2), sep = ' ')


#wprs
tt1 <-glm(maxIndurationTumour~wprs_alan + age_at_radiotherapy_start_yrs + factor(smoker) + factor(chemo) + factor(history_of_heart_disease) + BMI + b3radio_breast_ct_volume_cm3 + doseBED +  factor(b3radio_boost) + factor(b3post_operative_infection) + factor(diabetes) + factor(ra), data = STAT_prs_factors)
summary(tt1) 
ttt1 <- confint(tt1, level = 0.95)
a2 <- paste(format(round(ttt1[2,1], 3), nsmall = 3), format(round(ttt1[2,2], 3), nsmall = 3), sep=', ')
aa2<- paste('(',a2,')', sep = '')
aaa2 <- paste(format(round(as.numeric(summary(tt1)$coeff[2]), 3), nsmall = 3), aa2, format(round(summary(tt1)$coeff[2,4], 2), nsmall = 2), sep = ' ')

aaa1
aaa2

#prs percentile
tt1 <- glm(maxIndurationTumour~prs_precentile_alan + age_at_radiotherapy_start_yrs + factor(smoker) + factor(chemo) + factor(history_of_heart_disease) + BMI + b3radio_breast_ct_volume_cm3 + doseBED +  factor(b3radio_boost) + factor(b3post_operative_infection) + factor(diabetes) + factor(ra), data = STAT_prs_factors)
summary(tt1) 
ttt1 <- confint(tt1, level = 0.95)
a1 <- paste(format(round(ttt1[2,1], 3), nsmall = 3), format(round(ttt1[2,2], 3), nsmall = 3), sep=', ')
aa1<- paste('(',a1,')', sep = '')
aaa1 <- paste(format(round(as.numeric(summary(tt1)$coeff[2]), 3), nsmall = 3), aa1, format(round(summary(tt1)$coeff[2,4], 2), nsmall = 2), sep = ' ')


#wprs percentile
tt1 <- glm(maxIndurationTumour~wprs_precentile_alan + age_at_radiotherapy_start_yrs + factor(smoker) + factor(chemo) + factor(history_of_heart_disease) + BMI + b3radio_breast_ct_volume_cm3 + doseBED +  factor(b3radio_boost) + factor(b3post_operative_infection) + factor(diabetes) + factor(ra), data = STAT_prs_factors)
summary(tt1) 
ttt1 <- confint(tt1, level = 0.95)
a2 <- paste(format(round(ttt1[2,1], 3), nsmall = 3), format(round(ttt1[2,2], 3), nsmall = 3), sep=', ')
aa2<- paste('(',a2,')', sep = '')
aaa2 <- paste(format(round(as.numeric(summary(tt1)$coeff[2]), 3), nsmall = 3), aa2, format(round(summary(tt1)$coeff[2,4], 2), nsmall = 2), sep = ' ')

aaa1
aaa2



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



