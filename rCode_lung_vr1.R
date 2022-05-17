

############################################################################################
####
#### Script for calculating STAT from toxicity data
#### A McWilliam
#### vr1 7th May 2021
####
############################################################################################

### 1. Read data in and merge
########  Dataframes different sizes - toxicity entries multiple per patient
########  Will need date of treatment start/end to calculate time
### 2. select ID and colums of toxicit y to calculated STAT - timepoints?
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

### select ID and radiotherapy start data
radiotherapyStart <- Treat %>%
  select(SubjectId, radiotherapy_start_date)

### select out toxicity data of interest for STAT
toxicity <- Tox %>%
  select(SubjectId, date, cough, dyspnoea, broncho_pulmonary_haemorrhage, pneumonitis, pulmonary_fibrosis, esophagitis, dysphagia)

toxicity <- merge(toxicity, radiotherapyStart, by = 'SubjectId')
toxicity$monthStartTreat <- interval(ymd(as.Date(toxicity$radiotherapy_start_date)), ymd(as.Date(toxicity$date)))
toxicity$monthStartTreat = toxicity$monthStartTreat %/% months(1)

View(toxicity)

### need to set a time point for selecting highest toxicity score and filter
### keep highest recorded score for each toxicity 
months = 3
minMonths = 0

## pulmonary_fibrosis,  broncho_pulmonary_haemorrhage
## max_broncho_pulmonary_haemorrhage = max(broncho_pulmonary_haemorrhage), max_pulmonary_fibrosis = max(pulmonary_fibrosis),
toxicityFiltered <- toxicity %>%
  select(SubjectId, monthStartTreat, cough, dyspnoea, pneumonitis,  esophagitis, dysphagia) %>%
  filter(monthStartTreat <= months)  %>%
  filter(monthStartTreat > 1) %>%
  group_by(SubjectId) %>%
  summarise(maxMonths = max(monthStartTreat), max_cough = max(cough), max_dyspnoea = max(dyspnoea), max_pneumonitis = max(pneumonitis), max_esophagitis = max(esophagitis), max_dysphagia = max(dysphagia))

View(toxicityFiltered)

### check for complete cases and remove NA
toxicityFiltered$comp <- complete.cases(toxicityFiltered)
toxicityFiltered <- toxicityFiltered[toxicityFiltered$comp == TRUE,]
toxicityFiltered <- toxicityFiltered[-ncol(toxicityFiltered)]

view(dfSummary(toxicityFiltered))


### calculated the STAT for each patient
calcSTAT <- function(df){
  ### allocate storage for STAT values
  STATtmp <- matrix(NA, nrow = nrow(df), ncol = 1)
  
  ###calculate column mean and standard deviations
  col_Mean <- df %>% summarise_if(is.numeric, mean)
  col_SD <- df %>% summarise_if(is.numeric, sd)
  
  for(i in 1:nrow(df)){
    tmp = 0
    for(j in 1:length(df)){
      if(col_SD[j] != 0){
        tmp = tmp + (df[i,j] - col_Mean[j])/col_SD[j]
        }
    }
    tmp = tmp/length(df)
    STATtmp[i] = tmp
  }
  
  return(STATtmp)
}

STAT = calcSTAT(toxicityFiltered[ ,3:ncol(toxicityFiltered)])  ### removing the first two colums, studyID and months

### join STAT back on to data 
toxicityFilteredSTAT <- bind_cols(toxicityFiltered, t(as.data.frame(STAT)))
names(toxicityFilteredSTAT)[ncol(toxicityFilteredSTAT)] <- "STAT"
View(toxicityFilteredSTAT)

### select other clinical variables for inclusion in analysis


### save csv of data wkith STAT ready for analysis.
write.csv(toxicityFilteredSTAT, "C:/Users/alan_/Desktop/RAanalysis/REQUITEdata/processed/prostate.csv")


############################################################################################
### summary of patient characteristics
### summary STAT and plot histogram
###

summary(toxicityFilteredSTAT$STAT)

ggplot(data=toxicityFilteredSTAT, aes(STAT)) + 
  geom_histogram(breaks=seq(-0.7,2.5, by = 0.2),
                 col = "skyblue", fill = "lightblue") +
  labs(title = i, x = "STAT" ) +
  theme(panel.background = element_blank())

### save plots - change name
ggsave("C:/Users/alan_/Desktop/rheumotology/REQUITEdata/processed/figures/STATprostate.jpg")



############################################################################################
### association with polygenetic risk scores

load("C:/Users/alan_/Desktop/rheumotology/REQUITEdata/geneDose/data/REQUITE_PolygenicRiskScores_11May2021.RData")

View(requite_risk_scores)

### this has studyID
View(requite_risk_scores[[2]])
View(requite_risk_scores[["prs"]])
View(requite_risk_scores[["wprs"]])

prs <- cbind(requite_risk_scores[[2]], requite_risk_scores[["prs"]], requite_risk_scores[["wprs"]])
colnames(prs) <- c("SubjectId", "Site", "Country", "CancerType", "prs", "wprs")
View(prs)


STAT_prs <- merge(toxicityFilteredSTAT, prs, by = "SubjectId")
View(STAT_prs)

summary(STAT_prs$prs)
ggplot(data=STAT_prs, aes(prs)) + 
  geom_histogram(breaks=seq(50,95, by = 1),
                 col = "skyblue", fill = "lightblue") +
  labs(title = i, x = "prs" ) +
  theme(panel.background = element_blank())

summary(STAT_prs$wprs)
ggplot(data=STAT_prs, aes(wprs)) + 
  geom_histogram(breaks=seq(5,15, by = 0.2),
                 col = "skyblue", fill = "lightblue") +
  labs(title = i, x = "wprs" ) +
  theme(panel.background = element_blank())

plot(STAT_prs$STAT, STAT_prs$prs)
plot(STAT_prs$STAT, STAT_prs$wprs)


##########################################################################################
### need to select other patient factors

patFactors <- prosFactor %>%
  select(SubjectId, age_at_radiotherapy_start_yrs, smoker)

patTreat <- prosTreat %>%
  select( SubjectId, p3radio_number_fractions)

STAT_prs_factors <- merge(STAT_prs, patFactors, by = "SubjectId")
STAT_prs_factors <- merge(STAT_prs, patTreat, by = "SubjectId")
View(STAT_prs)

##########################################################################################

# + factor(smoker)

cor.test( STAT_prs$STAT,  STAT_prs$prs, method = c("spearman"))

t <- glm(STAT~prs + p3radio_number_fractions, data = STAT_prs_factors)
t <- glm(STAT~prs, data = STAT_prs)

ggplot(STAT_prs_factors, aes(x= STAT, y = wprs)) +
  geom_point() 

+
  geom_smooth()

plot(t, which = 1)
plot(t, which = 2)

t <- glm(STAT~wprs, data = STAT_prs)
summary(t)

summary(STAT_prs_factors$prs)

boxplot(STAT~p3radio_number_fractions, STAT_prs_factors)

tt <- glm(max_cough~prs, data = STAT_prs)
summary(tt)
tt <- glm(max_dyspnoea~prs, data = STAT_prs)
summary(tt)
boxplot(log(prs)~maxUrinary_incontinence, STAT_prs_factors)
tt <- glm(max_pneumonitis~prs, data = STAT_prs)
summary(tt)
tt <- glm(max_esophagitis~prs, data = STAT_prs)
summary(tt)
tt <- glm(max_dysphagia~prs, data = STAT_prs)
summary(tt)


#### is this useful for tables, check...
##tbl_regression(variable, conf.level = 0.95, exponentiate = TRUE)

