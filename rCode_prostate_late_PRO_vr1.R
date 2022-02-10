

############################################################################################
####
#### Script for calculating STAT from toxicity data using patient reported outcomes
#### Late toxicity analysis
#### A McWilliam
#### vr1 initial code for Brian's PRS 
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


############################################################################################
### Data cleaning and calculating the STAT


### read in prostate toxicities
### also need to know data of radiotherapy
prosTox <- read.csv("C:/Users/alan_/Desktop/rheumotology/REQUITEdata/Prostate/study/datasets/dataset5011.tsv", sep="\t", header=T )
View(prosTox)
prosTreat <- read.csv("C:/Users/alan_/Desktop/rheumotology/REQUITEdata/Prostate/study/datasets/dataset5009.tsv", sep="\t", header=T)
View(prosTreat)
prosFactor <- read.csv("C:/Users/alan_/Desktop/rheumotology/REQUITEdata/Prostate/study/datasets/dataset5010.tsv", sep="\t", header=T)
View(prosFactor)


prosPRO <- read.csv("C:/Users/alan_/Desktop/rheumotology/REQUITEdata/Prostate/study/datasets/dataset5003.tsv", sep="\t", header=T)
View(prosPRO)

### select ID and radiotherapy start data
radiotherapyStart <- prosTreat %>%
  select(SubjectId, p3radio_startdate)

### select out toxicity data of interest for STAT
### late toxicities of interest: proctities, rectal bleeding, urinary frequency, hematurary, retention/decreased stream
### do for patient and clinician reported

PRO <- prosPRO %>%
  select(SubjectId, event_date, P5a_q07_urinary_frequency, P5a_q14a_blood_urine, P5a_q08_urinary_flow, P5a_q14b_blood_bowels, P5a_q15_sticky_slimy_motions)


toxicityPRO <- merge(PRO, radiotherapyStart, by = 'SubjectId')
toxicityPRO$monthStartTreat <- interval(ymd(as.Date(toxicityPRO$p3radio_startdate)), ymd(as.Date(toxicityPRO$event_date)))
toxicityPRO$monthStartTreat = toxicityPRO$monthStartTreat %/% months(1)
View(toxicityPRO)


### need to set a time point for selecting highest toxicity score and filter
### keep highest recorded score for each toxicity 
minMonths = 0

######################################################################
###  PRO
######################################################################

## identify patients with baseline values
baselineCounts <- toxicityPRO %>%
  group_by(SubjectId) %>%
  filter(monthStartTreat == 0) %>%
  select(SubjectId)

## filter by baseline and subtract first entry
toxicityPROSubtract <- toxicityPRO %>%
  select(SubjectId, monthStartTreat, P5a_q07_urinary_frequency, P5a_q14a_blood_urine, P5a_q08_urinary_flow, P5a_q14b_blood_bowels, P5a_q15_sticky_slimy_motions) %>%
  filter(SubjectId %in% baselineCounts$SubjectId) %>%
  group_by(SubjectId) %>%
  mutate(P5a_q07_urinary_frequency_diff = P5a_q07_urinary_frequency - first(P5a_q07_urinary_frequency), P5a_q14a_blood_urine_diff = P5a_q14a_blood_urine - first(P5a_q14a_blood_urine), P5a_q08_urinary_flow_diff = P5a_q08_urinary_flow - first(P5a_q08_urinary_flow), P5a_q14b_blood_bowels_diff = P5a_q14b_blood_bowels - first(P5a_q14b_blood_bowels), P5a_q15_sticky_slimy_motions_diff = P5a_q15_sticky_slimy_motions - first(P5a_q15_sticky_slimy_motions)) %>%
  filter(monthStartTreat > minMonths) %>%
  summarise(maxMonths = max(monthStartTreat), maxUrinaryFrequency = max(P5a_q07_urinary_frequency_diff), maxHematuria = max(P5a_q14a_blood_urine_diff), maxDecreasedUrinaryStream = max(P5a_q08_urinary_flow_diff), maxRectalBleeding = max(P5a_q14b_blood_bowels_diff), maxRectalMucus = max(P5a_q15_sticky_slimy_motions_diff)) %>%
  select(SubjectId, maxMonths, maxUrinaryFrequency, maxHematuria, maxDecreasedUrinaryStream, maxRectalBleeding, maxRectalMucus)

toxicityPROSubtract[toxicityPROSubtract < 0] <- 0 
View(toxicityPROSubtract)

## patients without baseline
toxicityPRONoBaseline <- toxicityPRO %>%
  filter(!SubjectId %in% baselineCounts$SubjectId) %>%
  filter(monthStartTreat > minMonths) %>%
  group_by(SubjectId) %>%
  summarise(maxMonths = max(monthStartTreat), maxUrinaryFrequency = max(P5a_q07_urinary_frequency), maxHematuria = max(P5a_q14a_blood_urine), maxDecreasedUrinaryStream = max(P5a_q08_urinary_flow), maxRectalBleeding = max(P5a_q14b_blood_bowels), maxRectalMucus = max(P5a_q15_sticky_slimy_motions)) %>%
  select(SubjectId, maxMonths, maxUrinaryFrequency, maxHematuria, maxDecreasedUrinaryStream, maxRectalBleeding,  maxRectalMucus)
  
View(toxicityPRONoBaseline)

## join both together again
toxicityPROFiltered <- rbind(toxicityPROSubtract, toxicityPRONoBaseline)
View(toxicityPROFiltered)

###create temp data frame to display summary stats
toxicityPRO_summaryStats <- toxicityPROFiltered %>%
  select(maxUrinaryFrequency, maxHematuria, maxDecreasedUrinaryStream, maxRectalBleeding,  maxRectalMucus)
view(dfSummary(toxicityPRO_summaryStats))

## calculate STAT and join
t <- toxicityPROFiltered[ ,3:ncol(toxicityPROFiltered)]
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
toxicityPROFilteredSTAT <- bind_cols(toxicityPROFiltered, t(as.data.frame(STAT)))
names(toxicityPROFilteredSTAT)[ncol(toxicityPROFilteredSTAT)] <- "STAT"
View(toxicityPROFilteredSTAT)


### save csv of data wkith STAT ready for analysis.
##write.csv(toxicityFilteredSTAT, "C:/Users/alan_/Desktop/rheumotology/REQUITEdata/processed/prostate_acuteSTAT.csv")


############################################################################################
### summary of patient characteristics
### summary STAT and plot histogram
###

summary(toxicityPROFilteredSTAT$STAT)

ggplot(data=toxicityPROFilteredSTAT, aes(STAT)) + 
  geom_histogram(breaks=seq(-1,3, by = 0.1),
                 col = "skyblue", fill = "lightblue") +
  labs(title = i, x = "STAT" ) +
  theme(panel.background = element_blank())


### save plots - change name
ggsave("C:/Users/alan_/Desktop/rheumotology/REQUITEdata/processed/figures/STATprostateLatePRO.jpg")



############################################################################################
### association with polygenetic risk scores

load("C:/Users/alan_/Desktop/rheumotology/REQUITEdata/geneDose/data/REQUITE_PolygenicRiskScores_11May2021.RData")

View(requite_risk_scores)

### this has studyID
View(requite_risk_scores[[2]])
View(requite_risk_scores[["prs"]])
View(requite_risk_scores[["wprs"]])
View(requite_risk_scores[["gene_dose"]])

#tmp <- requite_risk_scores[["gene_dose"]][,1:1]
#tmp
#sum(tmp)

prs <- cbind(requite_risk_scores[[2]], requite_risk_scores[["prs"]], requite_risk_scores[["wprs"]])
colnames(prs) <- c("SubjectId", "Site", "Country", "CancerType", "prs", "wprs")
#View(prs)


PRO_STAT_prs <- merge(toxicityPROFilteredSTAT, prs, by = "SubjectId")
View(PRO_STAT_prs)

summary(PRO_STAT_prs$prs)
ggplot(data=PRO_STAT_prs, aes(prs)) + 
  geom_histogram(breaks=seq(50,95, by = 1),
                 col = "skyblue", fill = "lightblue") +
  labs(title = i, x = "prs" ) +
  theme(panel.background = element_blank())

summary(PRO_STAT_prs$wprs)
ggplot(data=PRO_STAT_prs, aes(wprs)) + 
  geom_histogram(breaks=seq(5,15, by = 0.2),
                 col = "skyblue", fill = "lightblue") +
  labs(title = i, x = "wprs" ) +
  theme(panel.background = element_blank())

plot(PRO_STAT_prs$STAT, PRO_STAT_prs$prs)
plot(PRO_STAT_prs$STAT, PRO_STAT_prs$wprs)


##########################################################################################
### need to select other patient factors
### adjusting for androgen-deprivation therapy, prior prostatectomy, age at treatment, and total BED
### + diabetes? 
### BED for acute tox - 3 also try 10 for normal tissue???


patFactors <- prosFactor %>%
  select(SubjectId, age_at_radiotherapy_start_yrs, smoker, diabetes, ra)

patTreat <- prosTreat %>%
  select( SubjectId, p3radio_number_fractions, p3radio_externalbeam_dose_Gy, p3radical_prostatectomy, p3hormone_therapy, p3hormone_therapy_length_months)


### need to calculate the BED prescribed
### BED = D x (1 + [d / (??/??)])
alphaBeta <- 10 #3
patTreat$doseFraction <- patTreat$p3radio_externalbeam_dose_Gy/patTreat$p3radio_number_fractions
patTreat$doseBED <-   patTreat$p3radio_externalbeam_dose_Gy  * (1 + (patTreat$doseFraction/alphaBeta))



PRO_STAT_prs_factors <- merge(PRO_STAT_prs, patFactors, by = "SubjectId")
PRO_STAT_prs_factors <- merge(PRO_STAT_prs_factors, patTreat, by = "SubjectId")

### clean the data
PRO_STAT_prs_factors <- PRO_STAT_prs_factors %>%
  filter(CancerType == 'Prostate')

View(PRO_STAT_prs_factors)

view(dfSummary(PRO_STAT_prs_factors))

##########################################################################################

#STAT_prs_factors$Country
#t2 <- glm(wprs~Country, data = STAT_prs_factors)
#summary(t2)


t <- glm(STAT~prs, data = PRO_STAT_prs_factors)
summary(t)
t <- glm(STAT~wprs, data = PRO_STAT_prs_factors)
summary(t)
confint(t, level = 0.95)


#### smoker 
### Country
t <- glm(STAT~wprs + age_at_radiotherapy_start_yrs + diabetes + p3radical_prostatectomy + p3hormone_therapy + doseBED + ra, data = PRO_STAT_prs_factors)
summary(t) 
t <- glm(STAT~prs + age_at_radiotherapy_start_yrs + diabetes + p3radical_prostatectomy + p3hormone_therapy + doseBED, data = STAT_prs_factors)
summary(t) 
AIC(t)

##t <- glm(log(STAT)~wprs + age_at_radiotherapy_start_yrs + diabetes +  p3radical_prostatectomy + p3hormone_therapy + doseBED, data = STAT_prs_factors)
###summary(t) 


resid(t)
plot(density(resid(t))) #A density plot
qqnorm(resid(t)) # A quantile normal plot - good for checking normality
qqline(resid(t))



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

t_residual <- glm(wprs~residual, data = STAT_residuals)
summary(t_residual)


#### residuals with individual varients
geneDose = requite_risk_scores[["gene_dose"]]
## get names of snps
snpNames<-rownames(geneDose)
snpNames


geneDose <- as.data.frame(t(as.matrix(geneDose)))
geneDose <- cbind(requite_risk_scores[[2]], geneDose)
View(geneDose)

## join geneDose with residuals
geneDose = rename(geneDose, SubjectId = Id)
STAT_residuals_geneDose <- merge(STAT_residuals, geneDose, by = "SubjectId")
View(STAT_residuals_geneDose)



model_stats<-matrix(ncol=6,nrow=length(snpNames))

for (i in 1:length(snpNames)) {
  formula<-as.formula( paste(snpNames[i], paste( "residual" ), sep=" ~ " ) )
  model<-glm(formula, data=STAT_residuals_geneDose)
  print(formula)
  
  model_stats[i,1]<-as.numeric(coef(model)[2])#beta
  model_stats[i,2:3]<-as.numeric(confint(model,"residual"))#upper-, lower-CI
  model_stats[i,4]<-(model_stats[i,3] - model_stats[i,2])/(2*1.96)#SE
  model_stats[i,5]<-1/(model_stats[i,4]^2)#weights
  model_stats[i,6]<-summary(model)$coeff[2,4]#p-value
  
  
  
}

colnames(model_stats)<-c("Beta","Upper CI","Lower CI","SE","Weights","P-Value")
rownames(model_stats)<-snpNames	
View(model_stats)
write.csv(model_stats, "C:\\Users\\alan_\\Desktop\\rheumotology\\prostateAcuteResultssnps.csv")
summary(model_stats)

###### associate prs and wprs with toxicity endpoints
toxEndPoints <- colnames(toxicityFiltered)
toxEndPoints <- toxEndPoints[-c(1, 2)]
toxEndPoints




model_stats_toxPRS<-matrix(ncol=6,nrow=(length(toxEndPoints)))

for (i in 1:length(toxEndPoints)) {
  tmp <- paste("factor(",toxEndPoints[i],")")
  formula<-as.formula( paste(tmp, paste( "prs" ), sep=" ~ " ) )
  print(formula)
  model<-glm(formula, family=binomial(link='logit'), data=STAT_residuals_geneDose)
  
  model_stats_toxPRS[i,1]<-as.numeric(coef(model)[2])#beta
  #print(model_stats_toxPRS[i,1])
  model_stats_toxPRS[i,2:3]<-as.numeric(confint(model,"prs"))#upper-, lower-CI
  #print(model_stats_toxPRS[i,2:3])
  model_stats_toxPRS[i,4]<-(model_stats[i,3] - model_stats[i,2])/(2*1.96)#SE
  model_stats_toxPRS[i,5]<-1/(model_stats[i,4]^2)#weights
  model_stats_toxPRS[i,6]<-summary(model)$coeff[2,4]#p-value
  
  
}

colnames(model_stats_toxPRS)<-c("Beta","Upper CI","Lower CI","SE","Weights","P-Value")
rownames(model_stats_toxPRS)<-toxEndPoints	
View(model_stats_toxPRS)
write.csv(model_stats_toxPRS, "C:\\Users\\alan_\\Desktop\\rheumotology\\prostateAcuteResults_toxoictyPRS.csv")
summary(model_stats_toxPRS)


######
model_stats_toxWPRS<-matrix(ncol=6,nrow=(length(toxEndPoints)))

for (i in 1:length(toxEndPoints)) {
  tmp <- paste("factor(",toxEndPoints[i],")")
  formula<-as.formula( paste(tmp, paste( "wprs" ), sep=" ~ " ) )
  print(formula)
  model<-glm(formula, family=binomial(link='logit'), data=STAT_residuals_geneDose)
  
  model_stats_toxWPRS[i,1]<-as.numeric(coef(model)[2])#beta
  #print(model_stats_toxPRS[i,1])
  model_stats_toxWPRS[i,2:3]<-as.numeric(confint(model,"wprs"))#upper-, lower-CI
  #print(model_stats_toxPRS[i,2:3])
  model_stats_toxWPRS[i,4]<-(model_stats[i,3] - model_stats[i,2])/(2*1.96)#SE
  model_stats_toxWPRS[i,5]<-1/(model_stats[i,4]^2)#weights
  model_stats_toxWPRS[i,6]<-summary(model)$coeff[2,4]#p-value
  
  
}

colnames(model_stats_toxWPRS)<-c("Beta","Upper CI","Lower CI","SE","Weights","P-Value")
rownames(model_stats_toxWPRS)<-toxEndPoints	
View(model_stats_toxWPRS)
write.csv(model_stats_toxWPRS, "C:\\Users\\alan_\\Desktop\\rheumotology\\prostateAcuteResults_toxoictyWPRS.csv")
summary(model_stats_toxWPRS)
 

####################################################

### look at patients with Rhematoid arthritis
summary(factor(STAT_prs_factors$ra))
tapply(STAT_prs_factors$wprs, STAT_prs_factors$ra, summary)
tapply(STAT_prs_factors$prs, STAT_prs_factors$ra, summary)

STAT_prs_factors$ra <- as.factor(STAT_prs_factors$ra)
ggplot(STAT_prs_factors, aes(x = ra, y = wprs)) + 
  geom_boxplot() + 
  theme_classic() +
  stat_compare_means(method = "wilcox.test", aes(group = ra, label = paste0("p = ",..p.format..)), label.x = 1.4, label.y = 6, size = 6)

ggplot(STAT_prs_factors, aes(x = ra, y = prs)) + 
  geom_boxplot() + 
  theme_classic() +
  stat_compare_means(method = "wilcox.test", aes(group = ra, label = paste0("p = ",..p.format..)), label.x = 1.4, label.y = 50, size = 6)


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
sarahList <- read.csv("C:\\Users\\alan_\\Desktop\\rheumotology\\listSarah\\REQUITE_prostate_STATacute_for_Alan.txt", sep = '\t')
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
