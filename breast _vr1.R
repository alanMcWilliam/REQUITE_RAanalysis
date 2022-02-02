
#### librarys needed
library(ggplot2)
library(dplyr)
library(summarytools)
library(lubridate)
library(ggpubr)

### read in breast data
### also need to know data of radiotherapy
breastTox <- read.csv("C:/Users/alan_/Desktop/rheumotology/REQUITEdata/Breast/study/datasets/dataset5016.tsv", sep="\t", header=T )
View(breastTox)
breastTreat <- read.csv("C:/Users/alan_/Desktop/rheumotology/REQUITEdata/Breast/study/datasets/dataset5018.tsv", sep="\t", header=T)
View(breastTreat)
breastFactor <- read.csv("C:/Users/alan_/Desktop/rheumotology/REQUITEdata/Breast/study/datasets/dataset5021.tsv", sep="\t", header=T)
View(breastFactor)

### select ID and radiotherapy start data
radiotherapyStart <- breastTreat %>%
  select(SubjectId, b3radio_breast_startdate)

### select out toxicity data of interest for STAT
### use combined STAT and individual

### breast atrophy, breast nipple retraction, breast edema, breast skin ulceration, breast telangiectasia (tumour bed),
### breast telangiectasia (outside tumour bed), breast skin induration (tumour bed), breast skin induration (outside tumour bed),
### Breast erythema, Breast arm lymphedema, breast skin hyperpigmentation, pneumonities, breast pain

### STAT: Telangiectasia, fibrosis (atrophyretraction), edema/induration

tox <- breastTox %>%
  select(SubjectId, date, breast_atrophy, breast_nipple_retraction, breast_oedema, breast_skin_ulceration, breast_telangiectasia_tumour_bed, breast_telangiectasia_outside_tumour_bed, breast_telangiectasia_outside_tumour_bed, breast_skin_induration_tumour_bed, breast_skin_induration_outside_tumour_bed, breast_erythema, breast_arm_lymphodema, breast_skin_hyperpigmentation, breast_pneumontis, breast_pain)

toxicityBreast <- merge(tox, radiotherapyStart, by = 'SubjectId')
toxicityBreast$monthStartTreat <- interval(ymd(as.Date(toxicityBreast$b3radio_breast_startdate)), ymd(as.Date(toxicityBreast$date)))
toxicityBreast$monthStartTreat = toxicityBreast$monthStartTreat %/% months(1)
View(toxicityBreast)

### need to set a time point for selecting highest toxicity score and filter
### keep highest recorded score for each toxicity 
minMonths = 0
months = 3
# filter(monthStartTreat <= months)  %>% ## can add below for acute toxicity

######################################################################
## identify patients with baseline values
baselineCounts <- toxicityBreast %>%
  group_by(SubjectId) %>%
  filter(monthStartTreat == 0) %>%
  select(SubjectId)

## filter by baseline and subtract first entry
toxicityBreastSubtract <- toxicityBreast %>%
  select(SubjectId, monthStartTreat, breast_atrophy, breast_nipple_retraction, breast_oedema, breast_skin_ulceration, breast_telangiectasia_tumour_bed, breast_telangiectasia_outside_tumour_bed, breast_skin_induration_tumour_bed, breast_skin_induration_outside_tumour_bed, breast_erythema, breast_arm_lymphodema, breast_skin_hyperpigmentation, breast_pneumontis, breast_pain) %>%
  filter(SubjectId %in% baselineCounts$SubjectId) %>%
  group_by(SubjectId) %>%
  mutate(breast_atrophy_diff = breast_atrophy - first(breast_atrophy), 
         breast_nipple_retraction_diff = breast_nipple_retraction - first(breast_nipple_retraction), 
         breast_oedema_diff = breast_oedema - first(breast_oedema), 
         breast_skin_ulceration_diff = breast_skin_ulceration - first(breast_skin_ulceration),
         breast_telangiectasia_tumour_bed_diff = breast_telangiectasia_tumour_bed - first(breast_telangiectasia_tumour_bed), 
         breast_telangiectasia_outside_tumour_bed_diff = breast_telangiectasia_outside_tumour_bed - first(breast_telangiectasia_outside_tumour_bed), 
         breast_skin_induration_tumour_bed_diff = breast_skin_induration_tumour_bed - first(breast_skin_induration_tumour_bed), 
         breast_skin_induration_outside_tumour_bed_diff = breast_skin_induration_outside_tumour_bed - first(breast_skin_induration_outside_tumour_bed), 
         breast_erythema_diff = breast_erythema - first(breast_erythema), 
         breast_arm_lymphodema_diff = breast_arm_lymphodema - first(breast_arm_lymphodema), 
         breast_skin_hyperpigmentation_diff = breast_skin_hyperpigmentation - first(breast_skin_hyperpigmentation),
         breast_pneumontis_diff = breast_pneumontis - first(breast_pneumontis),
         breast_pain_diff = breast_pain - first(breast_pain)) %>%
  filter(monthStartTreat > minMonths) %>%
  #filter(monthStartTreat <= months)  %>%
  summarise(maxMonths = max(monthStartTreat), 
            maxbreast_atrophy = max(breast_atrophy_diff), 
            maxbreast_nipple_retraction = max(breast_nipple_retraction_diff), 
            maxbreast_oedema = max(breast_oedema_diff), 
            maxbreast_skin_ulceration = max(breast_skin_ulceration_diff), 
            maxbreast_telangiectasia_tumour_bed = max(breast_telangiectasia_tumour_bed_diff),
            maxbreast_telangiectasia_outside_tumour_bed = max(breast_telangiectasia_outside_tumour_bed_diff), 
            maxbreast_skin_induration_tumour_bed = max(breast_skin_induration_tumour_bed_diff), 
            maxbreast_skin_induration_outside_tumour_bed = max(breast_skin_induration_outside_tumour_bed_diff), 
            maxbreast_erythema = max(breast_erythema_diff), 
            maxbreast_arm_lymphodema = max(breast_arm_lymphodema_diff), 
            maxbreast_skin_hyperpigmentation = max(breast_skin_hyperpigmentation_diff), 
            maxbreast_pneumontis = max(breast_pneumontis_diff), 
            maxbreast_pain = max(breast_pain_diff)) %>%
  select(SubjectId, maxMonths, maxbreast_atrophy, maxbreast_nipple_retraction, maxbreast_oedema, maxbreast_skin_ulceration, maxbreast_telangiectasia_tumour_bed, maxbreast_telangiectasia_outside_tumour_bed, maxbreast_skin_induration_tumour_bed, maxbreast_skin_induration_outside_tumour_bed, maxbreast_erythema, maxbreast_arm_lymphodema, maxbreast_skin_hyperpigmentation, maxbreast_pneumontis, maxbreast_pain)

toxicityBreastSubtract[toxicityBreastSubtract < 0] <- 0 
View(toxicityBreastSubtract)

## patients without baseline
toxicityBreastNoBaseline <- toxicityBreast %>%
  filter(!SubjectId %in% baselineCounts$SubjectId) %>%
  filter(monthStartTreat > minMonths) %>%
  #filter(monthStartTreat <= months)  %>%
  group_by(SubjectId) %>%
  summarise(maxMonths = max(monthStartTreat), 
            maxbreast_atrophy = max(breast_atrophy), 
            maxbreast_nipple_retraction = max(breast_nipple_retraction), 
            maxbreast_oedema = max(breast_oedema), 
            maxbreast_skin_ulceration = max(breast_skin_ulceration), 
            maxbreast_telangiectasia_tumour_bed = max(breast_telangiectasia_tumour_bed),
            maxbreast_telangiectasia_outside_tumour_bed = max(breast_telangiectasia_outside_tumour_bed), 
            maxbreast_skin_induration_tumour_bed = max(breast_skin_induration_tumour_bed), 
            maxbreast_skin_induration_outside_tumour_bed = max(breast_skin_induration_outside_tumour_bed), 
            maxbreast_erythema = max(breast_erythema), 
            maxbreast_arm_lymphodema = max(breast_arm_lymphodema), 
            maxbreast_skin_hyperpigmentation = max(breast_skin_hyperpigmentation), 
            maxbreast_pneumontis = max(breast_pneumontis), 
            maxbreast_pain = max(breast_pain)) %>%
  select(SubjectId, maxMonths, maxbreast_atrophy, maxbreast_nipple_retraction, maxbreast_oedema, maxbreast_skin_ulceration, maxbreast_telangiectasia_tumour_bed, maxbreast_telangiectasia_outside_tumour_bed, maxbreast_skin_induration_tumour_bed, maxbreast_skin_induration_outside_tumour_bed, maxbreast_erythema, maxbreast_arm_lymphodema, maxbreast_skin_hyperpigmentation, maxbreast_pneumontis, maxbreast_pain)

View(toxicityBreastNoBaseline)

## join both together again
toxicityBreastFiltered <- rbind(toxicityBreastSubtract, toxicityBreastNoBaseline)
View(toxicityBreastFiltered)

###create temp data frame to display summary stats
toxicityBreast_summaryStats <- toxicityBreastFiltered %>%
  select(SubjectId, maxMonths, maxbreast_atrophy, maxbreast_nipple_retraction, maxbreast_oedema, maxbreast_skin_ulceration, maxbreast_telangiectasia_tumour_bed, maxbreast_telangiectasia_outside_tumour_bed, maxbreast_skin_induration_tumour_bed, maxbreast_skin_induration_outside_tumour_bed, maxbreast_erythema, maxbreast_arm_lymphodema, maxbreast_skin_hyperpigmentation, maxbreast_pneumontis, maxbreast_pain)
stview(dfSummary(toxicityBreast_summaryStats))

## calculate STAT and join
t2 <- toxicityBreastFiltered[ ,3:ncol(toxicityBreastFiltered)]
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
toxicityBreastFilteredSTAT <- bind_cols(toxicityBreastFiltered, t(as.data.frame(STAT)))
names(toxicityBreastFilteredSTAT)[ncol(toxicityBreastFilteredSTAT)] <- "STAT"
View(toxicityBreastFilteredSTAT)

############################################################################################
### summary of patient characteristics
### summary STAT and plot histogram
###

summary(toxicityBreastFilteredSTAT$STAT)

ggplot(data=toxicityBreastFilteredSTAT, aes(STAT)) + 
  geom_histogram(breaks=seq(-1,3, by = 0.2),
                 col = "skyblue", fill = "lightblue") +
  labs(title = i, x = "STAT" ) +
  theme(panel.background = element_blank())

### save plots - change name
#ggsave("C:/Users/alan_/Desktop/rheumotology/REQUITEdata/processed/figures/STATprostateLatePRO.jpg")


############################################################################################
### association with polygenetic risk scores
### Brians scores
load("C:/Users/alan_/Desktop/rheumotology/REQUITEdata/geneDose/data/REQUITE_PolygenicRiskScores_11May2021.RData")

#View(requite_risk_scores)

### this has studyID
#View(requite_risk_scores[[2]])
#View(requite_risk_scores[["prs"]])
#View(requite_risk_scores[["wprs"]])
#View(requite_risk_scores[["gene_dose"]])

#tmp <- requite_risk_scores[["gene_dose"]][,1:1]
#tmp
#sum(tmp)

prs <- cbind(requite_risk_scores[[2]], requite_risk_scores[["prs"]], requite_risk_scores[["wprs"]])
colnames(prs) <- c("SubjectId", "Site", "Country", "CancerType", "prs", "wprs")
#View(prs)


Breast_STAT_prs <- merge(toxicityBreastFilteredSTAT, prs, by = "SubjectId")
View(Breast_STAT_prs)

summary(Breast_STAT_prs$prs)
ggplot(data=Breast_STAT_prs, aes(prs)) + 
  geom_histogram(breaks=seq(50,95, by = 1),
                 col = "skyblue", fill = "lightblue") +
  labs(title = i, x = "prs" ) +
  theme(panel.background = element_blank())

summary(Breast_STAT_prs$wprs)
ggplot(data=Breast_STAT_prs, aes(wprs)) + 
  geom_histogram(breaks=seq(5,15, by = 0.2),
                 col = "skyblue", fill = "lightblue") +
  labs(title = i, x = "wprs" ) +
  theme(panel.background = element_blank())


############################################################################################
### association with polygenetic risk scores
### Sarah's scores
sarahPRS <- read.csv("C:/Users/alan_/Desktop/rheumotology/listSarah/REQUITE_breast_lung_prostate_RA_PRS.txt", sep = "\t", header = F)
colnames(sarahPRS) <- c("SampleID", "prs_sarah", "wprs_sarah")
sarahPRS$SampleID <- gsub('.{7}$', '', sarahPRS$SampleID)
#View(sarahPRS)

sample <- read.csv("C:/Users/alan_/Desktop/rheumotology/REQUITEdata/Breast/study/datasets/dataset5039.tsv", sep="\t", header=T)
#View(sample)

sarahPRS <- merge(sarahPRS, sample, by = "SampleID")
#View(sarahPRS)

#compPRS <- merge(sarahPRS, prs, by = "SubjectId")
#View(compPRS)

### join for testing
breast_STAT_prs <- merge(Breast_STAT_prs, sarahPRS, by = "SubjectId")
#View(breast_STAT_prs)

summary(breast_STAT_prs$prs)
summary(breast_STAT_prs$prs_sarah)

summary(breast_STAT_prs$wprs)
summary(breast_STAT_prs$wprs_sarah)

ggplot(data = breast_STAT_prs) +
  geom_histogram(aes(x = prs_sarah), 
                 alpha=0.3, fill ="red",binwidth=2,position="dodge") +
  geom_histogram(aes(x = prs), 
                 alpha=0.3, fill ="green",binwidth=2,position="dodge") +
  labs(title = "", x = "prs") +
  theme(panel.background = element_blank())

ggplot(data = breast_STAT_prs) +
  geom_histogram(aes(x = wprs_sarah), 
                 alpha=0.3, fill ="red",binwidth=0.2,position="dodge") +
  geom_histogram(aes(x = wprs), 
                 alpha=0.3, fill ="green",binwidth=0.2,position="dodge") +
  labs(title = "", x = "wprs") +
  theme(panel.background = element_blank())

#######################################################


t <- glm(STAT~prs, data = breast_STAT_prs)
summary(t)

t <- glm(STAT~prs_sarah, data = breast_STAT_prs)
summary(t)

t <- glm(STAT~wprs, data = breast_STAT_prs)
summary(t)

t <- glm(STAT~wprs_sarah, data = breast_STAT_prs)
summary(t)



breast_STAT_prs$prs_precentile <- breast_STAT_prs$prs> quantile(breast_STAT_prs$prs, c(.90)) 
breast_STAT_prs$prs_precentile_sarah <- breast_STAT_prs$prs_sarah> quantile(breast_STAT_prs$prs_sarah, c(.90)) 

t <- glm(STAT~prs_precentile, data = breast_STAT_prs)
summary(t)
t <- glm(STAT~prs_precentile_sarah, data = breast_STAT_prs)
summary(t)

breast_STAT_prs$wprs_precentile <- breast_STAT_prs$wprs> quantile(breast_STAT_prs$wprs, c(.9)) 
breast_STAT_prs$wprs_precentile_sarah <- breast_STAT_prs$wprs_sarah> quantile(breast_STAT_prs$wprs_sarah, c(.9)) 

t <- glm(STAT~wprs_precentile, data = breast_STAT_prs)
summary(t)
t <- glm(STAT~wprs_precentile_sarah, data = breast_STAT_prs)
summary(t)

############################################################
###### associate prs and wprs with toxicity endpoints
toxEndPoints <- colnames(toxicityBreastFiltered)
toxEndPoints <- toxEndPoints[-c(1, 2)]
toxEndPoints


model_stats_breast<-matrix(ncol=6,nrow=(length(toxEndPoints)))

for (i in 1:length(toxEndPoints)) {
  tmp <- paste("factor(",toxEndPoints[i],")")
  formula<-as.formula( paste(tmp, paste( "prs_sarah" ), sep=" ~ " ) )
  print(formula)
  model<-glm(formula, family=binomial(link='logit'), data=breast_STAT_prs)
  
  model_stats_breast[i,1]<-as.numeric(coef(model)[2])#beta
  #print(model_stats_breast[i,1])
  model_stats_breast[i,2:3]<-as.numeric(confint(model,"prs_sarah"))#upper-, lower-CI
  #print(model_stats_breast[i,2:3])
  model_stats_breast[i,4]<-(model_stats_breast[i,3] - model_stats_breast[i,2])/(2*1.96)#SE
  model_stats_breast[i,5]<-1/(model_stats_breast[i,4]^2)#weights
  model_stats_breast[i,6]<-summary(model)$coeff[2,4]#p-value
}

colnames(model_stats_breast)<-c("Beta","Upper CI","Lower CI","SE","Weights","P-Value")
rownames(model_stats_breast)<-toxEndPoints	
View(model_stats_breast)
#write.csv(model_stats_breast, "C:\\Users\\alan_\\Desktop\\rheumotology\\prostateAcuteResults_toxoictyPRS.csv")
#summary(model_stats_breast)


######

model_stats_breast_w<-matrix(ncol=6,nrow=(length(toxEndPoints)))

for (i in 1:length(toxEndPoints)) {
  tmp <- paste("factor(",toxEndPoints[i],")")
  formula<-as.formula( paste(tmp, paste( "wprs_sarah" ), sep=" ~ " ) )
  print(formula)
  model<-glm(formula, family=binomial(link='logit'), data=breast_STAT_prs)
  
  model_stats_breast_w[i,1]<-as.numeric(coef(model)[2])#beta
  #print(model_stats_breast[i,1])
  model_stats_breast_w[i,2:3]<-as.numeric(confint(model,"wprs_sarah"))#upper-, lower-CI
  #print(model_stats_breast[i,2:3])
  model_stats_breast_w[i,4]<-(model_stats_breast_w[i,3] - model_stats_breast_w[i,2])/(2*1.96)#SE
  model_stats_breast_w[i,5]<-1/(model_stats_breast_w[i,4]^2)#weights
  model_stats_breast_w[i,6]<-summary(model)$coeff[2,4]#p-value
}

colnames(model_stats_breast_w)<-c("Beta","Upper CI","Lower CI","SE","Weights","P-Value")
rownames(model_stats_breast_w)<-toxEndPoints	
View(model_stats_breast_w)
#write.csv(model_stats_breast_w, "C:\\Users\\alan_\\Desktop\\rheumotology\\prostateAcuteResults_toxoictyPRS.csv")
#summary(model_stats_breast_w)



########################################################################

RA_tmp <- breastFactor %>%
  select(SubjectId, ra)

breast_STAT_prs_ra <- merge(breast_STAT_prs, RA_tmp, by = 'SubjectId')

summary(factor(breast_STAT_prs_ra$ra))
ggplot(breast_STAT_prs_ra, aes(x = factor(ra), y = wprs)) + 
  geom_boxplot() + 
  theme_classic() +
  stat_compare_means(method = "wilcox.test", aes(group = ra, label = paste0("p = ",..p.format..)), label.x = 1.4, label.y = 6, size = 6)

ggplot(breast_STAT_prs_ra, aes(x = factor(ra), y = wprs_sarah)) + 
  geom_boxplot() + 
  theme_classic() +
  stat_compare_means(method = "wilcox.test", aes(group = ra, label = paste0("p = ",..p.format..)), label.x = 1.4, label.y = 6, size = 6)


ggplot(breast_STAT_prs_ra, aes(x = factor(ra), y = prs)) + 
  geom_boxplot() + 
  theme_classic() +
  stat_compare_means(method = "wilcox.test", aes(group = ra, label = paste0("p = ",..p.format..)), label.x = 1.4, label.y = 50, size = 6)

ggplot(breast_STAT_prs_ra, aes(x = factor(ra), y = prs_sarah)) + 
  geom_boxplot() + 
  theme_classic() +
  stat_compare_means(method = "wilcox.test", aes(group = ra, label = paste0("p = ",..p.format..)), label.x = 1.4, label.y = 50, size = 6)

