

############################################################################################
####
#### A McWilliam
#### vr1 7th May 2021 first code calculating PRS
#### vr2 8th June 2021 updated to include multivarible analysis
#### 
#### code updated to incorperate PRS from Sarah and then my calculation of the PRS 
#### 
#### vr 3 updated for new PRS calculation
####
############################################################################################

### 1. Read data in and merge
########  Dataframes different sizes - toxicity entries multiple per patient
########  Will need date of treatment start/end to calculate time
### 2. select ID and colums of toxicity to calculated STAT - timepoints?
### 3. calculate STAT for each patient
### 4. select variables out for model building
### 5. analyse against gene dose - logistic regressions / multi-variable / other ideas...  bootstrapping data etc

### start with prostate



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
prosTox <- read.csv("C:/Users/alan_/Desktop/RAanalysis/REQUITEdata/Prostate/study/datasets/dataset5011.tsv", sep="\t", header=T )
View(prosTox)
prosTreat <- read.csv("C:/Users/alan_/Desktop/RAanalysis/REQUITEdata/Prostate/study/datasets/dataset5009.tsv", sep="\t", header=T)
View(prosTreat)
prosFactor <- read.csv("C:/Users/alan_/Desktop/RAanalysis/REQUITEdata/Prostate/study/datasets/dataset5010.tsv", sep="\t", header=T)
View(prosFactor)


prosPRO <- read.csv("C:/Users/alan_/Desktop/RAanalysis/REQUITEdata/Prostate/study/datasets/dataset5003.tsv", sep="\t", header=T)
View(prosPRO)

### select ID and radiotherapy start data
radiotherapyStart <- prosTreat %>%
  select(SubjectId, p3radio_startdate)


### select out toxicity data of interest for STAT
### Prostate: urinary_tract_obstruction, urinary_retention, proctitis, perforation, bowel_obstruction, bowel_fistula,bowel_stenosis, bowel_ulceration,flatus, fistula

PRO <- prosPRO %>%
  select(SubjectId, event_date, P5a_q07_urinary_frequency, P5a_q13_urine_night, P5a_q06_urgency_urine, P5a_q14a_blood_urine, P5a_q12_bladder_control, P5a_q08_urinary_flow, P5a_q03_pain_urine, P5a_q14b_med_blood_bowels, P5a_q11_bowel_control, P5a_q09_diarrhoea, P5a_q04_urgency_bowels, P5a_q05_urgency_bowels_unable, P5a_q02_pain_open_bowels, P5a_q15_sticky_slimy_motions, P5a_q10_constipation)


##toxicity <- prosTox %>%
##  select(SubjectId, date, urinary_frequency, urinary_urgency, haematuria, urinary_incontinence, diarrhoea, rectal_bleeding, management_of_sphincter_control)

toxicity <- merge(PRO, radiotherapyStart, by = 'SubjectId')
toxicity$monthStartTreat <- interval(ymd(as.Date(toxicity$p3radio_startdate)), ymd(as.Date(toxicity$event_date)))
toxicity$monthStartTreat = toxicity$monthStartTreat %/% months(1)

View(toxicity)

### 

### need to set a time point for selecting highest toxicity score and filter
### keep highest recorded score for each toxicity 
months = 50 #3 ## acute toxicity - 90 days
minMonths = 0

## identify patients with baseline values
baselineCounts <- toxicity %>%
  group_by(SubjectId) %>%
  filter(monthStartTreat == 0) %>%
  select(SubjectId)

## filter by baseline and subtract first entry
toxicitySubtract <- toxicity %>%
  select(SubjectId, monthStartTreat, P5a_q07_urinary_frequency, P5a_q13_urine_night, P5a_q06_urgency_urine, P5a_q14a_blood_urine, P5a_q12_bladder_control, P5a_q08_urinary_flow, P5a_q03_pain_urine, P5a_q14b_med_blood_bowels, P5a_q11_bowel_control, P5a_q09_diarrhoea, P5a_q04_urgency_bowels, P5a_q05_urgency_bowels_unable, P5a_q02_pain_open_bowels, P5a_q15_sticky_slimy_motions, P5a_q10_constipation) %>%
  filter(SubjectId %in% baselineCounts$SubjectId) %>%
  filter(monthStartTreat <= months)  %>%
  group_by(SubjectId) %>%
  mutate(P5a_q07_urinary_frequency_diff = P5a_q07_urinary_frequency - first(P5a_q07_urinary_frequency), P5a_q13_urine_night_diff = P5a_q13_urine_night - first(P5a_q13_urine_night), P5a_q06_urgency_urine_diff = P5a_q06_urgency_urine - first(P5a_q06_urgency_urine), P5a_q14a_blood_urine_diff = P5a_q14a_blood_urine - first(P5a_q14a_blood_urine), P5a_q12_bladder_control_diff = P5a_q12_bladder_control - first(P5a_q12_bladder_control), P5a_q08_urinary_flow_diff = P5a_q08_urinary_flow - first(P5a_q08_urinary_flow), P5a_q03_pain_urine_diff = P5a_q03_pain_urine - first(P5a_q03_pain_urine), P5a_q14b_med_blood_bowels_diff = P5a_q14b_med_blood_bowels - first(P5a_q14b_med_blood_bowels), P5a_q11_bowel_control_diff = P5a_q11_bowel_control - first(P5a_q11_bowel_control), P5a_q09_diarrhoea_diff = P5a_q09_diarrhoea - first(P5a_q09_diarrhoea), P5a_q04_urgency_bowels_diff = P5a_q04_urgency_bowels - first(P5a_q04_urgency_bowels), P5a_q05_urgency_bowels_unable_diff = P5a_q05_urgency_bowels_unable - first(P5a_q05_urgency_bowels_unable), P5a_q02_pain_open_bowels_diff = P5a_q02_pain_open_bowels - first(P5a_q02_pain_open_bowels), P5a_q15_sticky_slimy_motions_diff = P5a_q15_sticky_slimy_motions - first(P5a_q15_sticky_slimy_motions), P5a_q10_constipation_diff = P5a_q10_constipation - first(P5a_q10_constipation)) %>%
  filter(monthStartTreat > minMonths) %>%
  summarise(maxMonths = max(monthStartTreat), maxUrinaryFrequency = max(P5a_q07_urinary_frequency_diff), maxNocturia = max(P5a_q13_urine_night_diff), maxUrinaryUrgency = max(P5a_q06_urgency_urine_diff), maxHematuria = max(P5a_q14a_blood_urine_diff), maxUrinaryIncontenence = max(P5a_q12_bladder_control_diff), maxDecreasedUrinaryStream = max(P5a_q08_urinary_flow_diff), maxDysuria = max(P5a_q03_pain_urine_diff), maxRectalBleeding = max(P5a_q14b_med_blood_bowels_diff), maxGastrointenstinalIncontenence = max(P5a_q11_bowel_control_diff), maxDiarrhoea = max(P5a_q09_diarrhoea_diff), maxGastrointenstinalUrgency = max(P5a_q04_urgency_bowels_diff), maxTenesmus = max(P5a_q05_urgency_bowels_unable_diff), maxGastrointestinalPain = max(P5a_q02_pain_open_bowels_diff), maxRectalMucus = max(P5a_q15_sticky_slimy_motions_diff), maxConstipation = max(P5a_q10_constipation_diff)) %>%
  select(SubjectId, maxMonths, maxUrinaryFrequency, maxNocturia, maxUrinaryUrgency, maxHematuria, maxUrinaryIncontenence, maxDecreasedUrinaryStream, maxDysuria, maxGastrointenstinalIncontenence, maxDiarrhoea, maxGastrointenstinalUrgency, maxTenesmus, maxGastrointestinalPain, maxRectalMucus, maxConstipation)

toxicitySubtract[toxicitySubtract < 0] <- 0 
View(toxicitySubtract)


## patients without baseline
toxicityNoBaseline <- toxicity %>%
  filter(!SubjectId %in% baselineCounts$SubjectId) %>%
  filter(monthStartTreat <= months)  %>%
  filter(monthStartTreat > minMonths) %>%
  group_by(SubjectId) %>%
  summarise(maxMonths = max(monthStartTreat), maxUrinaryFrequency = max(P5a_q07_urinary_frequency), maxNocturia = max(P5a_q13_urine_night), maxUrinaryUrgency = max(P5a_q06_urgency_urine), maxHematuria = max(P5a_q14a_blood_urine), maxUrinaryIncontenence = max(P5a_q12_bladder_control), maxDecreasedUrinaryStream = max(P5a_q08_urinary_flow), maxDysuria = max(P5a_q03_pain_urine), maxRectalBleeding = max(P5a_q14b_med_blood_bowels), maxGastrointenstinalIncontenence = max(P5a_q11_bowel_control), maxDiarrhoea = max(P5a_q09_diarrhoea), maxGastrointenstinalUrgency = max(P5a_q04_urgency_bowels), maxTenesmus = max(P5a_q05_urgency_bowels_unable), maxGastrointestinalPain = max(P5a_q02_pain_open_bowels), maxRectalMucus = max(P5a_q15_sticky_slimy_motions), maxConstipation = max(P5a_q10_constipation)) %>%
  select(SubjectId, maxMonths, maxUrinaryFrequency, maxNocturia, maxUrinaryUrgency, maxHematuria, maxUrinaryIncontenence, maxDecreasedUrinaryStream, maxDysuria, maxGastrointenstinalIncontenence, maxDiarrhoea, maxGastrointenstinalUrgency, maxTenesmus, maxGastrointestinalPain, maxRectalMucus, maxConstipation)
  
View(toxicityNoBaseline)

## join both together again
toxicityFiltered <- rbind(toxicitySubtract, toxicityNoBaseline)
View(toxicityFiltered)


###create temp data frame to display summary stats
toxicity_summaryStats <- toxicityFiltered %>%
  select(maxUrinaryFrequency, maxNocturia, maxUrinaryUrgency, maxHematuria, maxUrinaryIncontenence, maxDecreasedUrinaryStream, maxDysuria, maxGastrointenstinalIncontenence, maxDiarrhoea, maxGastrointenstinalUrgency, maxTenesmus, maxGastrointestinalPain, maxRectalMucus, maxConstipation)

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
View(toxicityFilteredSTAT)

### select other clinical variables for inclusion in analysis


### save csv of data wkith STAT ready for analysis.
#write.csv(toxicityFilteredSTAT, "C:/Users/alan_/Desktop/RAanalysis/REQUITEdata/processed/prostate_acuteSTAT.csv")


############################################################################################
### summary of patient characteristics
### summary STAT and plot histogram
###

summary(toxicityFilteredSTAT$STAT)

ggplot(data=toxicityFilteredSTAT, aes(STAT)) + 
  geom_histogram(breaks=seq(-1,3, by = 0.1),
                 col = "skyblue", fill = "lightblue") +
  labs(title = i, x = "STAT" ) +
  theme(panel.background = element_blank())


ggplot(data=toxicityFilteredSTAT, aes(log10(STAT))) + 
  geom_histogram(breaks=seq(-3,3, by = 0.1),
                 col = "skyblue", fill = "lightblue") +
  labs(title = i, x = "STAT" ) +
  theme(panel.background = element_blank())


### save plots - change name
#ggsave("C:/Users/alan_/Desktop/RAanalysis/REQUITEdata/processed/figures/STATprostate.jpg")


##########################################################################################
### load in PRS + wPRS

alanPRS <- read.csv("C:/Users/alan_/Desktop/RAanalysis/calcPRS/PRS_all.csv", header = F)
alanPRS <- t(alanPRS)
alanPRS <- alanPRS[-1,]

colnames(alanPRS) <- c("SampleID", "prs_alan")
alanPRS <- as.data.frame(alanPRS)
alanPRS$SampleID <- as.numeric(alanPRS$SampleID)
alanPRS$prs_alan <- as.numeric(alanPRS$prs_alan)
#View(alanPRS)

alanWPRS <- read.csv("C:/Users/alan_/Desktop/RAanalysis/calcPRS/wPRS_all.csv", header = F)
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
sampleIDlink <- read.csv("C:/Users/alan_/Desktop/RAanalysis/REQUITEdata/Prostate/study/datasets/dataset5039.tsv", sep="\t", header=T)

##t <- merge(sampleIDlink, toxicityCROFilteredSTAT, by = "SubjectId")
STAT_prs <- merge(PRS_all, sampleIDlink, by = "SampleID")
STAT_prs <- merge(STAT_prs, toxicityFilteredSTAT, by = "SubjectId")
View(STAT_prs)



##########################################################################################
### need to select other patient factors
### adjusting for androgen-deprivation therapy, prior prostatectomy, age at treatment, and total BED
### + diabetes? 
### BED for acute tox - 3 also try 10 for normal tissue???


patFactors <- prosFactor %>%
  select(SubjectId, age_at_radiotherapy_start_yrs, smoker, diabetes, ra)

patTreat <- prosTreat %>%
  select( SubjectId, p3radio_number_fractions, p3radio_externalbeam_dose_Gy, p3radical_prostatectomy, p3hormone_therapy, p3hormone_therapy_length_months)

#tmp <- patTreat %>%
 # filter()

### need to calculate the BED prescribed
### BED = D x (1 + [d / (??/??)])
alphaBeta <- 10 #3
patTreat$doseFraction <- patTreat$p3radio_externalbeam_dose_Gy/patTreat$p3radio_number_fractions
patTreat$doseBED <-   patTreat$p3radio_externalbeam_dose_Gy  * (1 + (patTreat$doseFraction/alphaBeta))



STAT_prs_factors <- merge(STAT_prs, patFactors, by = "SubjectId")
STAT_prs_factors <- merge(STAT_prs_factors, patTreat, by = "SubjectId")

patCharct <- merge(patFactors, patTreat, by = "SubjectId")
view(dfSummary(patCharct))


### clean the data
STAT_prs_factors <- STAT_prs_factors %>%
  filter(CancerType == 'Prostate')

View(STAT_prs_factors)

view(dfSummary(STAT_prs_factors))


##########################################################################################
### analysing PRS and wPRS


t <- glm(STAT~prs_alan, data = STAT_prs_factors)
summary(t)
t <- glm(STAT~wprs_alan, data = STAT_prs_factors)
summary(t)

## + ra
t <- glm(STAT~prs_alan + age_at_radiotherapy_start_yrs + diabetes + p3radical_prostatectomy + p3hormone_therapy + doseBED, data = STAT_prs_factors)
summary(t) 
t <- glm(STAT~wprs_alan + age_at_radiotherapy_start_yrs + diabetes + p3radical_prostatectomy + p3hormone_therapy + doseBED, data = STAT_prs_factors)
summary(t) 


STAT_prs_factors$prs_precentile_alan <- STAT_prs_factors$prs_alan > quantile(STAT_prs_factors$prs_alan, c(.95)) 
t <- glm(STAT~prs_precentile_alan, data = STAT_prs_factors)
summary(t)

t <- glm(STAT~prs_precentile_alan + age_at_radiotherapy_start_yrs + diabetes + p3radical_prostatectomy + p3hormone_therapy + doseBED, data = STAT_prs_factors)
summary(t) 



STAT_prs_factors$wprs_precentile_alan <- STAT_prs_factors$wprs_alan > quantile(STAT_prs_factors$wprs_alan, c(.85)) 
t <- glm(STAT~wprs_precentile_alan, data = STAT_prs_factors)
summary(t)

t <- glm(STAT~wprs_precentile_alan + age_at_radiotherapy_start_yrs + diabetes + p3radical_prostatectomy + p3hormone_therapy + doseBED, data = STAT_prs_factors)
summary(t) 




##########################################################################################
##########################################################################################
##########################################################################################
### testing other aspects here 
### NEEDS UPDATING WITH NEW PRS + WPRS
STAT_prs_factors$Country
t2 <- glm(wprs~Country, data = STAT_prs_factors)
summary(t2)


t <- glm(STAT~prs, data = STAT_prs_factors)
summary(t)
t <- glm(STAT~wprs, data = STAT_prs_factors)
summary(t)
confint(t, level = 0.95)

t <- glm(log(STAT)~wprs, data = STAT_prs_factors)
summary(t)

#### smoker 
## country
t <- glm(STAT~wprs + age_at_radiotherapy_start_yrs + diabetes + p3radical_prostatectomy + p3hormone_therapy + doseBED + ra, data = STAT_prs_factors)
summary(t) 
tt <- exp(cbind(OR = coef(t), confint(t)))
tt

t <- glm(STAT~prs + age_at_radiotherapy_start_yrs + diabetes + p3radical_prostatectomy + p3hormone_therapy + doseBED, data = STAT_prs_factors)
summary(t) 
AIC(t)

##t <- glm(log(STAT)~wprs + age_at_radiotherapy_start_yrs + diabetes +  p3radical_prostatectomy + p3hormone_therapy + doseBED, data = STAT_prs_factors)
###summary(t) 


resid(t)
plot(density(resid(t))) #A density plot
qqnorm(resid(t)) # A quantile normal plot - good for checking normality
qqline(resid(t))


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

t2 <- glm(STAT~ age_at_radiotherapy_start_yrs + diabetes +  p3radical_prostatectomy + p3hormone_therapy + doseBED, data = STAT_prs_factors)
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
write.csv(model_stats, "C:\\Users\\alan_\\Desktop\\RAanalysis\\prostateResidualsSnp.csv")
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
