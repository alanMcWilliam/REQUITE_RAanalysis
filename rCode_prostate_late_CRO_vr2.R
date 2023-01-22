

############################################################################################
####
#### Script for calculating STAT from toxicity data using clinician reported outcomes
#### Late toxicity analysis
#### A McWilliam
#### vr1 initial code for Brian's PRS 
#### vr2 updated with new PRS calculation
#### 
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
### late toxicities of interest: proctities, rectal bleeding, urinary frequency, hematurary, retention/decreased stream
### do for patient and clinician reported

CRO <- prosTox %>%
  select(SubjectId, date, proctitis, rectal_bleeding, haematuria, urinary_frequency, urinary_retention)


toxicityCRO <- merge(CRO, radiotherapyStart, by = 'SubjectId')
toxicityCRO$monthStartTreat <- interval(ymd(as.Date(toxicityCRO$p3radio_startdate)), ymd(as.Date(toxicityCRO$date)))
toxicityCRO$monthStartTreat = toxicityCRO$monthStartTreat %/% months(1)
View(toxicityCRO)
### 

### need to set a time point for selecting highest toxicity score and filter
### keep highest recorded score for each toxicity 
minMonths = 3

######################################################################
###  CRO
######################################################################

## identify patients with baseline values
baselineCounts <- toxicityCRO %>%
  group_by(SubjectId) %>%
  filter(monthStartTreat == 0) %>%
  select(SubjectId)

## filter by baseline and subtract first entry
toxicityCROSubtract <- toxicityCRO %>%
  select(SubjectId, monthStartTreat, proctitis, rectal_bleeding, haematuria, urinary_frequency, urinary_retention) %>%
  filter(SubjectId %in% baselineCounts$SubjectId) %>%
  group_by(SubjectId) %>%
  mutate(proctitis_diff = proctitis - first(proctitis), rectal_bleeding_diff = rectal_bleeding - first(rectal_bleeding), haematuria_diff = haematuria - first(haematuria), urinary_frequency_diff = urinary_frequency - first(urinary_frequency), urinary_retention_diff = urinary_retention - first(urinary_retention)) %>%
  filter(monthStartTreat > minMonths) %>%
  summarise(maxMonths = max(monthStartTreat), maxProctitis = max(proctitis_diff), maxRectal_bleeding = max(rectal_bleeding_diff), maxHaematuria = max(haematuria_diff), maxUrinary_frequency = max(urinary_frequency_diff), maxUrinary_retention = max(urinary_retention_diff)) %>%
  select(SubjectId, maxMonths, maxProctitis, maxRectal_bleeding, maxHaematuria, maxUrinary_frequency, maxUrinary_retention)

toxicityCROSubtract[toxicityCROSubtract < 0] <- 0 
View(toxicityCROSubtract)

## patients without baseline
toxicityCRONoBaseline <- toxicityCRO %>%
  filter(!SubjectId %in% baselineCounts$SubjectId) %>%
  filter(monthStartTreat > minMonths) %>%
  group_by(SubjectId) %>%
  summarise(maxMonths = max(monthStartTreat), maxProctitis = max(proctitis), maxRectal_bleeding = max(rectal_bleeding), maxHaematuria = max(haematuria), maxUrinary_frequency = max(urinary_frequency), maxUrinary_retention = max(urinary_retention)) %>%
  select(SubjectId, maxMonths, maxProctitis, maxRectal_bleeding, maxHaematuria, maxUrinary_frequency, maxUrinary_retention)
View(toxicityCRONoBaseline)

## join both together again
toxicityCROFiltered <- rbind(toxicityCROSubtract, toxicityCRONoBaseline)
View(toxicityCROFiltered)

###create temp data frame to display summary stats
toxicityCRO_summaryStats <- toxicityCROFiltered %>%
  select(maxProctitis, maxRectal_bleeding, maxHaematuria, maxUrinary_frequency, maxUrinary_retention)
view(dfSummary(toxicityCRO_summaryStats))

## calculate STAT and join
t2 <- toxicityCROFiltered[ ,3:ncol(toxicityCROFiltered)]
col_Mean <- t2 %>% summarise_if(is.numeric, mean, na.rm = TRUE)
col_SD <- t2 %>% summarise_if(is.numeric, sd, na.rm = TRUE)
STAT <- matrix(NA, nrow = nrow(t2), ncol = 1)

for(i in 1:nrow(t2)){
  tmp2 = 0
  count2 = 0
  for(j in 1:length(t2)){
    if(!is.na(t2[i,j])){
      tmp2 = tmp2 + (t2[i,j] - col_Mean[j])/col_SD[j]
      count2 = count2 + 1
    }
  }
  tmp2 = tmp2/count2
  STAT[i] = tmp2
}
View(STAT)

### join STAT back on to data 
toxicityCROFilteredSTAT <- bind_cols(toxicityCROFiltered, t(as.data.frame(STAT)))
names(toxicityCROFilteredSTAT)[ncol(toxicityCROFilteredSTAT)] <- "STAT"
View(toxicityCROFilteredSTAT)


### save csv of data wkith STAT ready for analysis.
##write.csv(toxicityFilteredSTAT, "C:/Users/alan_/Desktop/rheumotology/REQUITEdata/processed/prostate_acuteSTAT.csv")


############################################################################################
### summary of patient characteristics
### summary STAT and plot histogram
###

summary(toxicityCROFilteredSTAT$STAT)

ggplot(data=toxicityCROFilteredSTAT, aes(STAT)) + 
  geom_histogram(breaks=seq(-1,4, by = 0.3),
                 col = "skyblue", fill = "lightblue") +
  labs(title = i, x = "STAT" ) +
  theme(panel.background = element_blank())


### save plots - change name
ggsave("C:/Users/alan_/Desktop/rheumotology/REQUITEdata/processed/figures/STATprostate.jpg")



##########################################################################################
### load in PRS + wPRS

alanPRS <- read.csv("C:/Users/alan_/Desktop/RAanalysis/calcPRS/PRS_all2.csv", header = F)
alanPRS <- t(alanPRS)
alanPRS <- alanPRS[-1,]

colnames(alanPRS) <- c("SampleID", "prs_alan")
alanPRS <- as.data.frame(alanPRS)
alanPRS$SampleID <- as.numeric(alanPRS$SampleID)
alanPRS$prs_alan <- as.numeric(alanPRS$prs_alan)
View(alanPRS)

alanWPRS <- read.csv("C:/Users/alan_/Desktop/RAanalysis/calcPRS/wPRS_all2.csv", header = F)
alanWPRS <- t(alanWPRS)
alanWPRS <- alanWPRS[-1,]

colnames(alanWPRS) <- c("SampleID", "wprs_alan")
alanWPRS <- as.data.frame(alanWPRS)
alanWPRS$SampleID <- as.numeric(alanWPRS$SampleID)
alanWPRS$wprs_alan <- as.numeric(alanWPRS$wprs_alan)
View(alanWPRS)

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
  labs(title = "", x = "prs") +
  theme(panel.background = element_blank())


## need to link back to patientID
sampleIDlink <- read.csv("C:/Users/alan_/Desktop/RAanalysis/REQUITEdata/Prostate/study/datasets/dataset5039.tsv", sep="\t", header=T)

CRO_STAT_prs <- merge(PRS_all, sampleIDlink, by = "SampleID")
CRO_STAT_prs <- merge(CRO_STAT_prs, toxicityCROFilteredSTAT, by = "SubjectId")
View(CRO_STAT_prs)


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



CRO_STAT_prs_factors <- merge(CRO_STAT_prs, patFactors, by = "SubjectId")
CRO_STAT_prs_factors <- merge(CRO_STAT_prs_factors, patTreat, by = "SubjectId")

### clean the data
#############
### this doesn't work with updated code - need to grab treatment site elsewhere to check
#############
CRO_STAT_prs_factors <- CRO_STAT_prs_factors %>%
  filter(CancerType == 'Prostate')

View(CRO_STAT_prs_factors)

view(dfSummary(CRO_STAT_prs_factors))

##########################################################################################


t <- glm(STAT~prs_alan, data = CRO_STAT_prs_factors)
summary(t)
t <- glm(STAT~wprs_alan, data = CRO_STAT_prs_factors)
summary(t)



t <- glm(STAT~wprs_alan + age_at_radiotherapy_start_yrs + diabetes + p3radical_prostatectomy + p3hormone_therapy + doseBED + ra, data = CRO_STAT_prs_factors)
summary(t) 
t <- glm(STAT~prs_alan + age_at_radiotherapy_start_yrs + diabetes + p3radical_prostatectomy + p3hormone_therapy + doseBED + ra, data = CRO_STAT_prs_factors)
summary(t) 

CRO_STAT_prs_factors$prs_precentile_alan <- CRO_STAT_prs_factors$prs_alan > quantile(CRO_STAT_prs_factors$prs_alan, c(.95)) 
t <- glm(STAT~prs_precentile_alan, data = CRO_STAT_prs_factors)
summary(t)
t <- glm(STAT~prs_precentile_alan + age_at_radiotherapy_start_yrs + diabetes + p3radical_prostatectomy + p3hormone_therapy + doseBED, data = CRO_STAT_prs_factors)
summary(t) 

CRO_STAT_prs_factors$wprs_precentile_alan <- CRO_STAT_prs_factors$wprs_alan > quantile(CRO_STAT_prs_factors$wprs_alan, c(.95)) 
t <- glm(STAT~wprs_precentile_alan, data = CRO_STAT_prs_factors)
summary(t)
t <- glm(STAT~wprs_precentile_alan + age_at_radiotherapy_start_yrs + diabetes + p3radical_prostatectomy + p3hormone_therapy + doseBED, data = CRO_STAT_prs_factors)
summary(t) 

##############################################################################################
##############################################################################################
##############################################################################################
##############################################################################################

### not using residual code here
##############################################################################################
##############################################################################################


##### work with residuals

t2 <- glm(STAT~ age_at_radiotherapy_start_yrs + diabetes +  p3radical_prostatectomy + p3hormone_therapy + doseBED, data = CRO_STAT_prs_factors)
summary(t2) 
AIC(t2)

resid(t2)
plot(density(resid(t2))) #A density plot
qqnorm(resid(t2)) # A quantile normal plot - good for checking normality
qqline(resid(t2))


residualsTmp <- t2$residuals 
residuals_present <- intersect(names(residualsTmp),rownames(CRO_STAT_prs_factors))
STAT_residuals <- CRO_STAT_prs_factors[residuals_present,]
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
genedose <- read.csv("C:/Users/alan_/Desktop/rheumotology/calcPRS/genedoses.csv", header = F)
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
write.csv(model_stats, "C:\\Users\\alan_\\Desktop\\rheumotology\\prostateResidualsSnp_lateCRO.csv")
summary(model_stats)



######################################################################
######################################################################
######################################################################
######################################################################
###### associate prs and wprs with toxicity endpoints
toxEndPoints <- colnames(toxicityCROFiltered)
toxEndPoints <- toxEndPoints[-c(1, 2)]
toxEndPoints


model_stats_toxPRS<-matrix(ncol=6,nrow=(length(toxEndPoints)))

for (i in 1:length(toxEndPoints)) {
  tmp <- paste("factor(",toxEndPoints[i],")")
  formula<-as.formula( paste(tmp, paste( "prs_alan" ), sep=" ~ " ) )
  print(formula)
  model<-glm(formula, family=binomial(link='logit'), data=CRO_STAT_prs_factors)
  
  model_stats_toxPRS[i,1]<-as.numeric(coef(model)[2])#beta
  #print(model_stats_toxPRS[i,1])
  model_stats_toxPRS[i,2:3]<-as.numeric(confint(model,"prs_alan"))#upper-, lower-CI
  #print(model_stats_toxPRS[i,2:3])
  model_stats_toxPRS[i,4]<-(model_stats_toxPRS[i,3] - model_stats_toxPRS[i,2])/(2*1.96)#SE
  model_stats_toxPRS[i,5]<-1/(model_stats_toxPRS[i,4]^2)#weights
  model_stats_toxPRS[i,6]<-summary(model)$coeff[2,4]#p-value
}

colnames(model_stats_toxPRS)<-c("Beta","Upper CI","Lower CI","SE","Weights","P-Value")
rownames(model_stats_toxPRS)<-toxEndPoints	
View(model_stats_toxPRS)
write.csv(model_stats_toxPRS, "C:\\Users\\alan_\\Desktop\\rheumotology\\CROprostateLateResults_PRS.csv")
summary(model_stats_toxPRS)


######
model_stats_toxWPRS<-matrix(ncol=6,nrow=(length(toxEndPoints)))

for (i in 1:length(toxEndPoints)) {
  tmp <- paste("factor(",toxEndPoints[i],")")
  formula<-as.formula( paste(tmp, paste( "wprs_alan" ), sep=" ~ " ) )
  print(formula)
  model<-glm(formula, family=binomial(link='logit'), data=CRO_STAT_prs_factors)
  
  model_stats_toxWPRS[i,1]<-as.numeric(coef(model)[2])#beta
  #print(model_stats_toxPRS[i,1])
  model_stats_toxWPRS[i,2:3]<-as.numeric(confint(model,"wprs_alan"))#upper-, lower-CI
  #print(model_stats_toxPRS[i,2:3])
  model_stats_toxWPRS[i,4]<-(model_stats_toxWPRS[i,3] - model_stats_toxWPRS[i,2])/(2*1.96)#SE
  model_stats_toxWPRS[i,5]<-1/(model_stats_toxWPRS[i,4]^2)#weights
  model_stats_toxWPRS[i,6]<-summary(model)$coeff[2,4]#p-value
}

colnames(model_stats_toxWPRS)<-c("Beta","Upper CI","Lower CI","SE","Weights","P-Value")
rownames(model_stats_toxWPRS)<-toxEndPoints	
View(model_stats_toxWPRS)
write.csv(model_stats_toxWPRS, "C:\\Users\\alan_\\Desktop\\rheumotology\\prostateAcuteResults_CROprostateLateResults_WPRS.csv")
summary(model_stats_toxWPRS)


## percentiles

model_stats_toxPRS<-matrix(ncol=6,nrow=(length(toxEndPoints)))

for (i in 1:length(toxEndPoints)) {
  tmp <- paste("factor(",toxEndPoints[i],")")
  formula<-as.formula( paste(tmp, paste( "prs_precentile_alan" ), sep=" ~ " ) )
  print(formula)
  model<-glm(formula, family=binomial(link='logit'), data=CRO_STAT_prs_factors)
  
  model_stats_toxPRS[i,1]<-as.numeric(coef(model)[2])#beta
  #print(model_stats_toxPRS[i,1])
  #model_stats_toxPRS[i,2:3]<-as.numeric(confint(model,"prs_precentile_alan"))#upper-, lower-CI
  #print(model_stats_toxPRS[i,2:3])
  model_stats_toxPRS[i,4]<-(model_stats_toxPRS[i,3] - model_stats_toxPRS[i,2])/(2*1.96)#SE
  model_stats_toxPRS[i,5]<-1/(model_stats_toxPRS[i,4]^2)#weights
  model_stats_toxPRS[i,6]<-summary(model)$coeff[2,4]#p-value
}

colnames(model_stats_toxPRS)<-c("Beta","Upper CI","Lower CI","SE","Weights","P-Value")
rownames(model_stats_toxPRS)<-toxEndPoints	
View(model_stats_toxPRS)
write.csv(model_stats_toxPRS, "C:\\Users\\alan_\\Desktop\\rheumotology\\CROprostateLateResults_PRSpercentile.csv")
summary(model_stats_toxPRS)


######
model_stats_toxWPRS<-matrix(ncol=6,nrow=(length(toxEndPoints)))

for (i in 1:length(toxEndPoints)) {
  tmp <- paste("factor(",toxEndPoints[i],")")
  formula<-as.formula( paste(tmp, paste( "wprs_precentile_alan" ), sep=" ~ " ) )
  print(formula)
  model<-glm(formula, family=binomial(link='logit'), data=CRO_STAT_prs_factors)
  
  model_stats_toxWPRS[i,1]<-as.numeric(coef(model)[2])#beta
  #print(model_stats_toxPRS[i,1])
  #model_stats_toxWPRS[i,2:3]<-as.numeric(confint(model,"wprs_precentile_alan"))#upper-, lower-CI
  #print(model_stats_toxPRS[i,2:3])
  model_stats_toxWPRS[i,4]<-(model_stats_toxWPRS[i,3] - model_stats_toxWPRS[i,2])/(2*1.96)#SE
  model_stats_toxWPRS[i,5]<-1/(model_stats_toxWPRS[i,4]^2)#weights
  model_stats_toxWPRS[i,6]<-summary(model)$coeff[2,4]#p-value
}

colnames(model_stats_toxWPRS)<-c("Beta","Upper CI","Lower CI","SE","Weights","P-Value")
rownames(model_stats_toxWPRS)<-toxEndPoints	
View(model_stats_toxWPRS)
write.csv(model_stats_toxWPRS, "C:\\Users\\alan_\\Desktop\\rheumotology\\prostateAcuteResults_CROprostateLateResults_WPRSpercentile.csv")
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
