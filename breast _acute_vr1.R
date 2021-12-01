
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

######################
### here
######################
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
