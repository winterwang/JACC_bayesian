#' ---
#' title: "JACC study Milk intake and stroke mortality analysis"
#' author: "Chaochen Wang"
#' date: "2019-12-20 created, `r Sys.Date()` updated"
#' output:
#'   pdf_document: 
#'     toc: true
#'     toc_depth: 3  # upto three depths of headings (specified by #, ## and ###)
#'     number_sections: true
#'     keep_tex: true
#'     latex_engine: xelatex 
#' header-includes: 
#'   - \usepackage{bookmark} 
#'   - \usepackage{xltxtra} 
#'   - \usepackage{zxjatype} 
#'   - \usepackage[ipaex]{zxjafont} 
#' ---
#' # Read in the data

library(readr)
library(tidyverse)
library(lubridate) # for dealing with date time data 

MILK <- read_csv("../data/StrokeMilk.csv", 
                     progress = show_progress(), 
                     col_types = cols(.default = "c"))

MILK %>% 
  filter(tr_age > 39 & tr_age < 80) %>% 
  group_by(tr_sex) %>% 
  summarise(n= n()) %>%
  mutate(rel.freq = paste0(round(100 * n/sum(n), 2), "%"))  

#' # delete subjects outside of age range ------------------------------------ 

MILK_0 <- MILK %>%
  filter(tr_age > 39 & tr_age < 80)

#' # define total stroke mortality --------------------------------


MILK_0 <- MILK_0 %>% 
  mutate(Tot_Stroke = if_else(grepl("I6[0-9][0-9]|I6[0-9]",  
                                    ICD10), "I60_9", 
                        if_else(!is.na(ICD10), "other_death", 
                                      "Alive/Censor"))) 
MILK_0%>% 
  group_by(tr_sex, Tot_Stroke) %>%
  summarise(n= n()) %>%
  mutate(rel.freq = paste0(round(100 * n/sum(n), 2), "%"))

#' # define different type of stroke mortality/CVD ?---------------------------

#' I60  Nontraumatic subarachnoid hemorrhage
#' 
#' I61  Nontraumatic intracerebral hemorrhage
#' 
#' I62  Other and unspecified nontraumatic intracranial hemorrhage
#' 
#' I63  Cerebral infarction
#' 
#' I65  Occlusion and stenosis of precerebral arteries, not resulting in cerebral infarction
#' 
#' I66  Occlusion and stenosis of cerebral arteries, not resulting in cerebral infarction
#' 
#' I67  Other cerebrovascular diseases
#' 
#' I68  Cerebrovascular disorders in diseases classified elsewhere
#' 
#' I69  Sequelae of cerebrovascular disease

MILK_0 <- MILK_0 %>% 
  mutate(HemoStroke = if_else(grepl("I6[0-2][0-9]|I6[0-2]",  
                                    ICD10), "I60_2", 
                              if_else(!is.na(ICD10), "other_death", 
                                      "Alive/Censor"))) %>% 
  mutate(IscheStroke = if_else(grepl("I63[0-9]|I63",  
                                    ICD10), "I63", 
                              if_else(!is.na(ICD10), "other_death", 
                                      "Alive/Censor"))) %>%
  mutate(CHD = if_else(grepl("I2[0-5][0-9]|I2[0-5]",  
                                     ICD10), "I20_5", 
                               if_else(!is.na(ICD10), "other_death", 
                                       "Alive/Censor"))) %>%
  mutate(HeartF = if_else(grepl("I50[0-9]|I50",  
                                     ICD10), "I50", 
                               if_else(!is.na(ICD10), "other_death", 
                                       "Alive/Censor")))


MILK_0%>% 
  group_by(tr_sex, HemoStroke) %>%
  summarise(n= n()) %>%
  mutate(rel.freq = paste0(round(100 * n/sum(n), 2), "%"))
MILK_0%>% 
  group_by(tr_sex, IscheStroke) %>%
  summarise(n= n()) %>%
  mutate(rel.freq = paste0(round(100 * n/sum(n), 2), "%"))
MILK_0%>% 
  group_by(tr_sex, CHD) %>%
  summarise(n= n()) %>%
  mutate(rel.freq = paste0(round(100 * n/sum(n), 2), "%"))
MILK_0%>% 
  group_by(tr_sex, HeartF) %>%
  summarise(n= n()) %>%
  mutate(rel.freq = paste0(round(100 * n/sum(n), 2), "%"))

#' # Define milk intake ---------------------------------------------
#' 

MILK_0 <- MILK_0 %>% 
  mutate(Milk_fre = as.numeric(MILK)) %>% 
  mutate(Milk_fre = as.factor(Milk_fre)) %>% 
  mutate(Mlkfre = fct_collapse(Milk_fre,
                               Never = "1",
                               Mon1_2 = "2", 
                               Wek1_2 = "3",
                               Wek3_4 = "4",
                               Daily  = "5")) %>% 
  mutate(MlkLogi = fct_collapse(Mlkfre,
                                Never = "Never", 
                                Drinker = c("Mon1_2", "Wek1_2", "Wek3_4", "Daily")))
MILK_0 %>% 
  group_by(tr_sex, Mlkfre) %>%
  summarise(n= n()) %>%
  mutate(rel.freq = paste0(round(100 * n/sum(n), 2), "%"))
MILK_0 %>% 
  group_by(tr_sex, MlkLogi) %>%
  summarise(n= n()) %>%
  mutate(rel.freq = paste0(round(100 * n/sum(n), 2), "%"))

#' # Calculate person-years

MILK_0 <- MILK_0 %>% 
  mutate(Age = as.numeric(tr_age)) %>% 
  mutate(Agegrp = cut(as.numeric(tr_age), c(30, 45, 55, 65, 75, 80), right = FALSE)) %>% 
  mutate(followpy = as.numeric(actual)/365.25) 

#' # Identify potential confounders: smoking, alcohol intake, BMI, DM/HYT/MI/APO/Cancer history, Exercise, Energy intake, Sleep duration, vegetable/fru/gretea/cofe intake, school education
#' 
#' 
MILK_0 <-  MILK_0 %>% 
  mutate(Smoking = replace_na(SM1, "unknown")) %>% 
  mutate(Smoking = as_factor(Smoking)) %>% 
  mutate(Smoking = fct_recode(Smoking, Never = "3", Past = "2", Current = "1")) %>% 
  mutate(Smoking = factor(Smoking, levels = c("Never", "Past", "Current", "unknown"))) %>%  # Smoking
  mutate(Alc_Fre = if_else(as.numeric(DR1) >= 2, "Never or past", 
                           if_else(as.numeric(DR1F) == 1, "Daily", 
                                   if_else(as.numeric(DR1F) == 4, "< 1/week", 
                                           if_else((as.numeric(DR1F) == 2) | (as.numeric(DR1F) == 3), 
                                                   "1-4 /week", "Unknown"))))) %>% 
  mutate(Alc_Fre = fct_explicit_na(Alc_Fre, na_level = "unknown")) %>% 
  mutate(BMI = as.numeric(wt10)/(as.numeric(ht10)^2) * 100000) %>% # define BMI groups
  mutate(BMIgrp = cut(BMI, breaks = c(14, 18.5, 25, 30, 40), right = FALSE)) %>% 
  mutate(BMIgrp = as.character(BMIgrp)) %>% 
  replace_na(list(BMIgrp = "unknown")) %>% 
  mutate(BMIgrp = factor(BMIgrp, levels = c("[18.5,25)",
                                            "[14,18.5)",
                                            "[25,30)", 
                                            "[30,40)", "unknown"))) %>%
  mutate(DM_hist = if_else(as.numeric(p_DM) > 1, TRUE, FALSE)) %>% 
  replace_na(list(DM_hist = "unknown")) %>% # recode DM history status
  mutate(HT_hist = if_else(as.numeric(p_HT) > 1, TRUE, FALSE)) %>% 
  replace_na(list(HT_hist = "unknown")) %>% # recode hyt history status
  mutate(MI_hist = if_else(as.numeric(p_MI) > 1, TRUE, FALSE)) %>% 
  replace_na(list(MI_hist = "unknown")) %>% # recode MI history status
  mutate(APO_hist = if_else(as.numeric(p_APO) > 1, TRUE, FALSE)) %>% 
  replace_na(list(APO_hist = "unknown")) %>% # recode APO history status
  mutate(KID_hist = if_else(as.numeric(p_KID) > 1, TRUE, FALSE)) %>% 
  replace_na(list(KID_hist = "unknown")) %>% # recode KID history status
  mutate(LIV_hist = if_else(as.numeric(p_APO) > 1, TRUE, FALSE)) %>% 
  replace_na(list(LIV_hist = "unknown")) %>% # recode LIV history status
  mutate(Can_hist = if_else(as.numeric(p_can1) > 1 | 
                              as.numeric(p_can2) > 1, TRUE, FALSE)) %>% 
  replace_na(list(Can_hist = "unknown")) %>% # recode LIV history status
  mutate(Exercise = as.numeric(sport) != 4) %>% # define exercise habits
  mutate(Exercise = as.character(Exercise)) %>% 
  replace_na(list(Exercise = "unknown")) %>% 
  mutate(Exercise = factor(Exercise, levels = c("FALSE", "TRUE", "unknown"))) %>% 
  mutate(Exercise = fct_recode(Exercise, 
                               "> 1h/w" = "TRUE", 
                               "Almost0" = "FALSE", 
                               unknown   = "unknown")) %>% 
  mutate(Engy = log(as.numeric(ENERGY))) %>% 
  mutate(Sleep = as.numeric(SLEEP)/10) %>% 
  mutate(Slepgrp = cut(Sleep, breaks = c(0, 6.9, 7.9, 8.9, 23), right = FALSE)) %>% 
  mutate(Slepgrp = as.character(Slepgrp)) %>% 
  replace_na(list(Slepgrp = "unknown")) %>% 
  mutate(Slepgrp = factor(Slepgrp, levels = c("[0,6.9)",
                                            "[6.9,7.9)",
                                            "[7.9,8.9)", 
                                            "[8.9,23)", "unknown"))) %>% 
  mutate(Spi = as.factor(SPI)) %>% # define vegetable intake
  mutate(Spi = fct_collapse(Spi, 
                            unknown = "X", 
                            daily = "5",
                            Thre4tw = "4",
                            One2tw = "3",
                            Less1tm = c("1", "2"))) %>% 
  mutate(Spi = fct_explicit_na(Spi, na_level = "unknown")) %>% 
  mutate(Fru = as.factor(FRU)) %>% # define fruit intake
  mutate(Fru = fct_collapse(Fru, 
                            unknown = "X",
                            daily = "5",
                            Thre4tw = "4",
                            One2tw = "3",
                            Less1tm = c("1", "2"))) %>% 
  mutate(Fru = fct_explicit_na(Fru, na_level = "unknown")) %>% 
  mutate(Gretea = as.factor(GreTEA1)) %>% # define greentea intake
  mutate(Gretea = fct_collapse(Gretea, 
                               unknown = "X", 
                               Thre3tw = "2", 
                               Thre3tw = "3",
                               Thre3tw = "4",
                               Never = "5", 
                               daily = "1")) %>% 
  mutate(Gretea = fct_explicit_na(Gretea, na_level = "unknown"))  %>% 
  mutate(Cofe = as.factor(COFE)) %>% # define greentea intake
  mutate(Cofe = fct_collapse(Cofe, 
                               unknown = "X", 
                               Thre3tw = "2", 
                               Thre3tw = "3",
                               Thre3tw = "4",
                               Never = "5", 
                               daily = "1")) %>% 
  mutate(Cofe = fct_explicit_na(Cofe, na_level = "unknown"))  %>% 
  mutate(Educ = as.numeric(MILK_0$SCHOOL)) %>% 
  mutate(Educgrp = cut(Educ, breaks = c(0, 18, 70), right = FALSE)) %>% 
  mutate(Educgrp = as.character(Educgrp)) %>% 
  replace_na(list(Educgrp = "unknown")) %>% 
  mutate(Educgrp = factor(Educgrp, levels = c("[0,18)",
                                              "[18,70)",
                                              "unknown"))) %>% # Define menopause for women
  mutate(Menopause = if_else(!is.na(MENO_AGE)& tr_sex == "2", TRUE,  # define menopause
                             if_else(as.numeric(tr_age) >= 50 & tr_sex == "2", 
                                     TRUE, FALSE)))

  
MILK_0 %>% 
  group_by(tr_sex, Smoking) %>% 
  summarise (n= n()) %>%
  mutate(rel.freq = paste0(round(100 * n/sum(n), 2), "%"))  %>% 
  print(n=Inf)
MILK_0 %>% 
  group_by(tr_sex, Alc_Fre) %>% 
  summarise (n= n()) %>%
  mutate(rel.freq = paste0(round(100 * n/sum(n), 2), "%"))  %>% 
  print(n=Inf)
MILK_0 %>% 
  group_by(tr_sex, BMIgrp) %>% 
  summarise (n= n()) %>%
  mutate(rel.freq = paste0(round(100 * n/sum(n), 2), "%"))  %>% 
  print(n=Inf)
MILK_0 %>% 
  group_by(tr_sex, DM_hist) %>% 
  summarise (n= n()) %>%
  mutate(rel.freq = paste0(round(100 * n/sum(n), 2), "%"))  %>% 
  print(n=Inf)
MILK_0 %>% 
  group_by(tr_sex, HT_hist) %>% 
  summarise (n= n()) %>%
  mutate(rel.freq = paste0(round(100 * n/sum(n), 2), "%"))  %>% 
  print(n=Inf)
MILK_0 %>% 
  group_by(tr_sex, MI_hist) %>% 
  summarise (n= n()) %>%
  mutate(rel.freq = paste0(round(100 * n/sum(n), 2), "%"))  %>% 
  print(n=Inf)
MILK_0 %>% 
  group_by(tr_sex, APO_hist) %>% 
  summarise (n= n()) %>%
  mutate(rel.freq = paste0(round(100 * n/sum(n), 2), "%"))  %>% 
  print(n=Inf)
MILK_0 %>% 
  group_by(tr_sex, KID_hist) %>% 
  summarise (n= n()) %>%
  mutate(rel.freq = paste0(round(100 * n/sum(n), 2), "%"))  %>% 
  print(n=Inf)
MILK_0 %>% 
  group_by(tr_sex, LIV_hist) %>% 
  summarise (n= n()) %>%
  mutate(rel.freq = paste0(round(100 * n/sum(n), 2), "%"))  %>% 
  print(n=Inf)
MILK_0 %>% 
  group_by(tr_sex, Can_hist) %>% 
  summarise (n= n()) %>%
  mutate(rel.freq = paste0(round(100 * n/sum(n), 2), "%"))  %>% 
  print(n=Inf)
MILK_0 %>% 
  group_by(tr_sex, Exercise) %>% 
  summarise (n= n()) %>%
  mutate(rel.freq = paste0(round(100 * n/sum(n), 2), "%"))  %>% 
  print(n=Inf)
MILK_0 %>% 
  group_by(tr_sex, Slepgrp) %>% 
  summarise (n= n()) %>%
  mutate(rel.freq = paste0(round(100 * n/sum(n), 2), "%"))  %>% 
  print(n=Inf)
MILK_0 %>% 
  group_by(tr_sex, Spi) %>% 
  summarise (n= n()) %>%
  mutate(rel.freq = paste0(round(100 * n/sum(n), 2), "%"))  %>% 
  print(n=Inf)
MILK_0 %>% 
  group_by(tr_sex, Fru) %>% 
  summarise (n= n()) %>%
  mutate(rel.freq = paste0(round(100 * n/sum(n), 2), "%"))  %>% 
  print(n=Inf)
MILK_0 %>% 
  group_by(tr_sex, Gretea) %>% 
  summarise (n= n()) %>%
  mutate(rel.freq = paste0(round(100 * n/sum(n), 2), "%"))  %>% 
  print(n=Inf)
MILK_0 %>% 
  group_by(tr_sex, Cofe) %>% 
  summarise (n= n()) %>%
  mutate(rel.freq = paste0(round(100 * n/sum(n), 2), "%"))  %>% 
  print(n=Inf)
MILK_0 %>% 
  group_by(tr_sex, Educgrp) %>% 
  summarise (n= n()) %>%
  mutate(rel.freq = paste0(round(100 * n/sum(n), 2), "%"))  %>% 
  print(n=Inf)
MILK_0 %>% 
  group_by(tr_sex, Menopause) %>% 
  summarise (n= n()) %>%
  mutate(rel.freq = paste0(round(100 * n/sum(n), 2), "%"))  %>% 
  print(n=Inf)


# 02-04 AREA 地区(施設番号+地区番号)
# – touhoku: (1, 2, 3, 4, 17, 29)
# – kanto: (5, 6, 8, 9, 11, 13, 31)
# – chubu: (15, 18)
# – kinki: (10, 20, 21, 22, 24)
# – chugoku: (25, 26)
# – kyushiu: (27, 30)

MILK_0 <- MILK_0 %>% 
  mutate(areano = as.numeric(areano)) %>% 
  mutate(Area = if_else(areano %in% c(11, 22, 23, 24, 41, 30,
                                      170, 178, 179, 298, 299), "Touhoku", 
                  if_else(areano %in% c(51, 61, 81, 91, 92, 93, 
                                          110, 130, 311), "Kanto", 
                    if_else(areano %in% c(151, 181), "Chubu", 
                      if_else(areano %in% c(100, 108, 109, 201, 211, 212, 213, 
                                214, 221, 241, 242, 243), "Kinki", 
                        if_else(areano %in% c(250, 261), "Chugoku", 
                          if_else(areano %in% c(271, 272, 273, 274, 300, 301, 302, 303, 304, 
                            305, 306, 307, 308, 309), "Kyushiu", "else"))))))) %>% 
  mutate(Area = factor(Area))


#' # Exclusion: history of stroke, cancer, MI, angina pectoris, other ischemic heart disease (ICD9)
#' 410-414  Ischemic Heart Disease
#' 
#' 415-417  Diseases Of Pulmonary Circulation
#' 
#' 420-429  Other Forms Of Heart Disease


MILK_0 %>% 
  group_by(tr_sex, APO_hist) %>% 
  summarise (n= n()) %>%
  mutate(rel.freq = paste0(round(100 * n/sum(n), 2), "%"))  %>% 
  print(n=Inf)


MILK_0 %>% 
  group_by(tr_sex, Can_hist) %>% 
  summarise (n= n()) %>%
  mutate(rel.freq = paste0(round(100 * n/sum(n), 2), "%"))  %>% 
  print(n=Inf)

MILK_0 %>% 
  group_by(tr_sex, MI_hist) %>% 
  summarise (n= n()) %>%
  mutate(rel.freq = paste0(round(100 * n/sum(n), 2), "%"))  %>% 
  print(n=Inf)


MILK_0 <- MILK_0 %>% 
  mutate(p_Oth1 = as.numeric(p_oth1c)) %>% 
  mutate(p_Oth2 = as.numeric(p_oth2c)) %>% 
  mutate(IscheHeart = if_else((p_Oth1 >=410 & p_Oth1 <=414) | 
                                (p_Oth2 >=410 & p_Oth2 <=414), TRUE, FALSE)) %>% 
  replace_na(list(IscheHeart = "unknown")) %>% # recode IscheHeart history status
  mutate(OtheHeart = if_else((p_Oth1 >=420 & p_Oth1 <=429) | 
                                (p_Oth2 >=420 & p_Oth2 <=429), TRUE, FALSE)) %>% 
  replace_na(list(OtheHeart = "unknown")) #%>% # recode Otherheart history status
  

MILK_0 %>% 
  group_by(tr_sex, IscheHeart) %>% 
  summarise (n= n()) %>%
  mutate(rel.freq = paste0(round(100 * n/sum(n), 2), "%"))  %>% 
  print(n=Inf)
MILK_0 %>% 
  group_by(tr_sex, OtheHeart) %>% 
  summarise (n= n()) %>%
  mutate(rel.freq = paste0(round(100 * n/sum(n), 2), "%"))  %>% 
  print(n=Inf)

MData <- MILK_0 %>%
  filter(APO_hist != "TRUE" & IscheHeart != "TRUE" & 
           OtheHeart != "TRUE" & Can_hist != "TRUE" & MI_hist != "TRUE" & !is.na(Mlkfre)) %>% 
  select(Area, Age, Agegrp, tr_sex, Tot_Stroke, HemoStroke, IscheStroke, CHD, HeartF, MlkLogi, 
         Mlkfre, followpy, Smoking, Alc_Fre, BMI, BMIgrp, DM_hist, HT_hist, KID_hist, 
         LIV_hist, Exercise, Engy, ENERGY, Sleep, Slepgrp, Spi, Fru, Gretea, Cofe, Educ,
         Educgrp, Menopause)

# data preparation done

MData_men <- MData %>% 
  filter(tr_sex == "1")
MData_fem <- MData %>% 
  filter(tr_sex == "2")

#' ## before entering the analyses ordered, we need to explore by preliminary analyses
#' 
#' 
# Number of subjects, number of cases, person years 
# by frequency

MData_men %>%
  group_by(Mlkfre) %>%
  summarise(pyear = sum(followpy), n = n()) %>% 
  mutate_if(is.numeric, format, 2)
MData_fem %>%
  group_by(Mlkfre) %>%
  summarise(pyear = sum(followpy), n = n()) %>% 
  mutate_if(is.numeric, format, 2)

epiDisplay::tabpct(MData_men$Mlkfre, MData_men$Tot_Stroke, 
                   percent = "row", graph = FALSE)
epiDisplay::tabpct(MData_fem$Mlkfre, MData_fem$Tot_Stroke, 
                   percent = "row", graph = FALSE)


###################################################################################################
## survival object
###################################################################################################
library(survival)
library(ggplot2)
library(survminer)
library(cowplot)
library(ggsci)

# in Men
su_obj_men <- Surv(MData_men$followpy, MData_men$Tot_Stroke == "I60_9")
# in Women
su_obj_fem <- Surv(MData_fem$followpy, MData_fem$Tot_Stroke == "I60_9")


###################################################################################################
## Kaplan-Meier plots and log rank tests for TotStroke and Milk intake  (frequency)
###################################################################################################


#+ fig1, echo=FALSE, eval=TRUE, fig.height=6.5, fig.width=7, message=FALSE, warning=FALSE, fig.cap="Kaplan-Meier survival curves for total stroke mortality by drinking frequency (P value was obtained from log-rank tests) in Men.",fig.align='center',  out.width='100%', cache=TRUE

surv.Drfreq <- with(MData_men, survfit(su_obj_men ~ Mlkfre))

ggsurv_Drfreq <- ggsurvplot(surv.Drfreq, censor = F, xlab = "Time (years)", conf.int = T, 
                            conf.int.style = "step",  # customize style of confidence intervals
                            surv.median.line = "none", ylab = "Survival function",
                            legend.labs = c("Never or past", "< 1/week", "1-2 /week", "3-4 /week", "Daily"), 
                            ggtheme = theme_bw(), palette = "npg",
                            pval = TRUE, data = MData_men, # pval.method = TRUE
                            risk.table = T,
                            risk.table.y.text.col = T, # colour risk table text annotations.
                            risk.table.y.text = FALSE,  pval.coord = c(0, 0.972),
                            ylim = c(0.93, 1.0), xlim = c(0, 22)) 
# show bars instead of names in text annotations)  

ggsurv_Drfreq <- ggpar(
  ggsurv_Drfreq,
  font.title    = c(14, "bold", "black"),         
  # font.subtitle = c(15, "bold.italic", "purple"), 
  #  font.caption  = c(14, "plain", "orange"),        
  font.x        = c(14, "bold", "black"),          
  font.y        = c(14, "bold", "black"),      
  font.xtickslab = c(13, "bold", "black"),
  #  font.ytickslab = c(13, "bold", "black"),
  legend = "bottom", 
  font.legend  = c(13, "bold", "black"),
  legend.title = ""
)
ggsurv_Drfreq

#+ fig2, echo=FALSE, eval=TRUE, fig.height=6.5, fig.width=7, message=FALSE, warning=FALSE, fig.cap="Kaplan-Meier survival curves for total stroke mortality by drinking frequency (P value was obtained from log-rank tests) in Women.",fig.align='center',  out.width='100%', cache=TRUE

surv.Drfreq <- with(MData_fem, survfit(su_obj_fem ~ Mlkfre))

ggsurv_Drfreq <- ggsurvplot(surv.Drfreq, censor = F, xlab = "Time (years)", conf.int = T, 
                            conf.int.style = "step",  # customize style of confidence intervals
                            surv.median.line = "none", ylab = "Survival function",
                            legend.labs = c("Never or past", "< 1/week", "1-2 /week", "3-4 /week", "Daily"), 
                            ggtheme = theme_bw(), palette = "npg",
                            pval = TRUE, data = MData_fem, # pval.method = TRUE
                            risk.table = T,
                            risk.table.y.text.col = T, # colour risk table text annotations.
                            risk.table.y.text = FALSE,  pval.coord = c(0, 0.972),
                            ylim = c(0.93, 1.0), xlim = c(0, 22)) 
# show bars instead of names in text annotations)  

ggsurv_Drfreq <- ggpar(
  ggsurv_Drfreq,
  font.title    = c(14, "bold", "black"),         
  # font.subtitle = c(15, "bold.italic", "purple"), 
  #  font.caption  = c(14, "plain", "orange"),        
  font.x        = c(14, "bold", "black"),          
  font.y        = c(14, "bold", "black"),      
  font.xtickslab = c(13, "bold", "black"),
  #  font.ytickslab = c(13, "bold", "black"),
  legend = "bottom", 
  font.legend  = c(13, "bold", "black"),
  legend.title = ""
)
ggsurv_Drfreq

#' # In Men
#' ## Model0

SurvM0 <-  coxph(su_obj_men ~ Mlkfre, 
                 data = MData_men)

library("broom")
tidy(SurvM0, exponentiate = TRUE, conf.int = TRUE) %>% 
  knitr::kable(.)



#' ## Model1 

SurvM1 <-  coxph(su_obj_men ~ Mlkfre + Age + strata(Agegrp), 
                 data = MData_men)

tidy(SurvM1, exponentiate = TRUE, conf.int = TRUE) %>% 
  knitr::kable(.)

#' ## Model2 

SurvM2 <-  coxph(su_obj_men ~ Mlkfre + Age + strata(Agegrp) + Smoking + Alc_Fre + 
                   BMIgrp + DM_hist + HT_hist + KID_hist + LIV_hist + Exercise + 
                   Slepgrp + Spi + Fru + Cofe + Educgrp + Gretea, 
                 data = MData_men)

tidy(SurvM2, exponentiate = TRUE, conf.int = TRUE) %>% 
  knitr::kable(.)

#' # In women 
#' ## Model0

SurvM0 <-  coxph(su_obj_fem ~ Mlkfre, 
                 data = MData_fem)

library("broom")
tidy(SurvM0, exponentiate = TRUE, conf.int = TRUE) %>% 
  knitr::kable(.)



#' ## Model1 

SurvM1 <-  coxph(su_obj_fem ~ Mlkfre + Age + strata(Agegrp), 
                 data = MData_fem)

tidy(SurvM1, exponentiate = TRUE, conf.int = TRUE) %>% 
  knitr::kable(.)

#' ## Model2 

SurvM2 <-  coxph(su_obj_fem ~ Mlkfre + Age + strata(Agegrp) + Smoking + Alc_Fre + 
                   BMIgrp + DM_hist + HT_hist + KID_hist + LIV_hist + Exercise + 
                   Slepgrp + Spi + Fru + Cofe + Educgrp + Gretea + Menopause, 
                 data = MData_fem)

tidy(SurvM2, exponentiate = TRUE, conf.int = TRUE) %>% 
  knitr::kable(.)

#' # Cause specific: HemoStroke

# in Men
su_obj_men <- Surv(MData_men$followpy, MData_men$HemoStroke == "I60_2")
# in Women
su_obj_fem <- Surv(MData_fem$followpy, MData_fem$HemoStroke == "I60_2")


#' # In Men
#' ## Model0

SurvM0 <-  coxph(su_obj_men ~ Mlkfre, 
                 data = MData_men)

library("broom")
tidy(SurvM0, exponentiate = TRUE, conf.int = TRUE) %>% 
  knitr::kable(.)



#' ## Model1 

SurvM1 <-  coxph(su_obj_men ~ Mlkfre + Age + strata(Agegrp), 
                 data = MData_men)

tidy(SurvM1, exponentiate = TRUE, conf.int = TRUE) %>% 
  knitr::kable(.)

#' ## Model2 

SurvM2 <-  coxph(su_obj_men ~ Mlkfre + Age + strata(Agegrp) + Smoking + Alc_Fre + 
                   BMIgrp + DM_hist + HT_hist + KID_hist + LIV_hist + Exercise + 
                   Slepgrp + Spi + Fru + Cofe + Educgrp + Gretea, 
                 data = MData_men)

tidy(SurvM2, exponentiate = TRUE, conf.int = TRUE) %>% 
  knitr::kable(.)

#' # In women 
#' ## Model0

SurvM0 <-  coxph(su_obj_fem ~ Mlkfre, 
                 data = MData_fem)

library("broom")
tidy(SurvM0, exponentiate = TRUE, conf.int = TRUE) %>% 
  knitr::kable(.)



#' ## Model1 

SurvM1 <-  coxph(su_obj_fem ~ Mlkfre + Age + strata(Agegrp), 
                 data = MData_fem)

tidy(SurvM1, exponentiate = TRUE, conf.int = TRUE) %>% 
  knitr::kable(.)

#' ## Model2 

SurvM2 <-  coxph(su_obj_fem ~ Mlkfre + Age + strata(Agegrp) + Smoking + Alc_Fre + 
                   BMIgrp + DM_hist + HT_hist + KID_hist + LIV_hist + Exercise + 
                   Slepgrp + Spi + Fru + Cofe + Educgrp + Gretea + Menopause, 
                 data = MData_fem)

tidy(SurvM2, exponentiate = TRUE, conf.int = TRUE) %>% 
  knitr::kable(.)



#' # Cause specific: IscheStroke

# in Men
su_obj_men <- Surv(MData_men$followpy, MData_men$IscheStroke == "I63")
# in Women
su_obj_fem <- Surv(MData_fem$followpy, MData_fem$IscheStroke == "I63")


#' # In Men
#' ## Model0

SurvM0 <-  coxph(su_obj_men ~ Mlkfre, 
                 data = MData_men)

library("broom")
tidy(SurvM0, exponentiate = TRUE, conf.int = TRUE) %>% 
  knitr::kable(.)



#' ## Model1 

SurvM1 <-  coxph(su_obj_men ~ Mlkfre + Age + strata(Agegrp), 
                 data = MData_men)

tidy(SurvM1, exponentiate = TRUE, conf.int = TRUE) %>% 
  knitr::kable(.)

#' ## Model2 

SurvM2 <-  coxph(su_obj_men ~ Mlkfre + Age + strata(Agegrp) + Smoking + Alc_Fre + 
                   BMIgrp + DM_hist + HT_hist + KID_hist + LIV_hist + Exercise + 
                   Slepgrp + Spi + Fru + Cofe + Educgrp + Gretea, 
                 data = MData_men)

tidy(SurvM2, exponentiate = TRUE, conf.int = TRUE) %>% 
  knitr::kable(.)

#' # In women 
#' ## Model0

SurvM0 <-  coxph(su_obj_fem ~ Mlkfre, 
                 data = MData_fem)

library("broom")
tidy(SurvM0, exponentiate = TRUE, conf.int = TRUE) %>% 
  knitr::kable(.)



#' ## Model1 

SurvM1 <-  coxph(su_obj_fem ~ Mlkfre + Age + strata(Agegrp), 
                 data = MData_fem)

tidy(SurvM1, exponentiate = TRUE, conf.int = TRUE) %>% 
  knitr::kable(.)

#' ## Model2 

SurvM2 <-  coxph(su_obj_fem ~ Mlkfre + Age + strata(Agegrp) + Smoking + Alc_Fre + 
                   BMIgrp + DM_hist + HT_hist + KID_hist + LIV_hist + Exercise + 
                   Slepgrp + Spi + Fru + Cofe + Educgrp + Gretea + Menopause, 
                 data = MData_fem)

tidy(SurvM2, exponentiate = TRUE, conf.int = TRUE) %>% 
  knitr::kable(.)

#' # Cause specific: CHD

# in Men
su_obj_men <- Surv(MData_men$followpy, MData_men$CHD == "I20_5")
# in Women
su_obj_fem <- Surv(MData_fem$followpy, MData_fem$CHD == "I20_5")


#' # In Men
#' ## Model0

SurvM0 <-  coxph(su_obj_men ~ Mlkfre, 
                 data = MData_men)

library("broom")
tidy(SurvM0, exponentiate = TRUE, conf.int = TRUE) %>% 
  knitr::kable(.)



#' ## Model1 

SurvM1 <-  coxph(su_obj_men ~ Mlkfre + Age + strata(Agegrp), 
                 data = MData_men)

tidy(SurvM1, exponentiate = TRUE, conf.int = TRUE) %>% 
  knitr::kable(.)

#' ## Model2 

SurvM2 <-  coxph(su_obj_men ~ Mlkfre + Age + strata(Agegrp) + Smoking + Alc_Fre + 
                   BMIgrp + DM_hist + HT_hist + KID_hist + LIV_hist + Exercise + 
                   Slepgrp + Spi + Fru + Cofe + Educgrp + Gretea, 
                 data = MData_men)

tidy(SurvM2, exponentiate = TRUE, conf.int = TRUE) %>% 
  knitr::kable(.)

#' # In women 
#' ## Model0

SurvM0 <-  coxph(su_obj_fem ~ Mlkfre, 
                 data = MData_fem)

library("broom")
tidy(SurvM0, exponentiate = TRUE, conf.int = TRUE) %>% 
  knitr::kable(.)



#' ## Model1 

SurvM1 <-  coxph(su_obj_fem ~ Mlkfre + Age + strata(Agegrp), 
                 data = MData_fem)

tidy(SurvM1, exponentiate = TRUE, conf.int = TRUE) %>% 
  knitr::kable(.)

#' ## Model2 

SurvM2 <-  coxph(su_obj_fem ~ Mlkfre + Age + strata(Agegrp) + Smoking + Alc_Fre + 
                   BMIgrp + DM_hist + HT_hist + KID_hist + LIV_hist + Exercise + 
                   Slepgrp + Spi + Fru + Cofe + Educgrp + Gretea + Menopause, 
                 data = MData_fem)

tidy(SurvM2, exponentiate = TRUE, conf.int = TRUE) %>% 
  knitr::kable(.)


#' # Cause specific: HeartF

# in Men
su_obj_men <- Surv(MData_men$followpy, MData_men$HeartF == "I50")
# in Women
su_obj_fem <- Surv(MData_fem$followpy, MData_fem$HeartF == "I50")


#' # In Men
#' ## Model0

SurvM0 <-  coxph(su_obj_men ~ Mlkfre, 
                 data = MData_men)

library("broom")
tidy(SurvM0, exponentiate = TRUE, conf.int = TRUE) %>% 
  knitr::kable(.)



#' ## Model1 

SurvM1 <-  coxph(su_obj_men ~ Mlkfre + Age + strata(Agegrp), 
                 data = MData_men)

tidy(SurvM1, exponentiate = TRUE, conf.int = TRUE) %>% 
  knitr::kable(.)

#' ## Model2 

SurvM2 <-  coxph(su_obj_men ~ Mlkfre + Age + strata(Agegrp) + Smoking + Alc_Fre + 
                   BMIgrp + DM_hist + HT_hist + KID_hist + LIV_hist + Exercise + 
                   Slepgrp + Spi + Fru + Cofe + Educgrp + Gretea, 
                 data = MData_men)

tidy(SurvM2, exponentiate = TRUE, conf.int = TRUE) %>% 
  knitr::kable(.)

#' # In women 
#' ## Model0

SurvM0 <-  coxph(su_obj_fem ~ Mlkfre, 
                 data = MData_fem)

library("broom")
tidy(SurvM0, exponentiate = TRUE, conf.int = TRUE) %>% 
  knitr::kable(.)



#' ## Model1 

SurvM1 <-  coxph(su_obj_fem ~ Mlkfre + Age + strata(Agegrp), 
                 data = MData_fem)

tidy(SurvM1, exponentiate = TRUE, conf.int = TRUE) %>% 
  knitr::kable(.)

#' ## Model2 

SurvM2 <-  coxph(su_obj_fem ~ Mlkfre + Age + strata(Agegrp) + Smoking + Alc_Fre + 
                   BMIgrp + DM_hist + HT_hist + KID_hist + LIV_hist + Exercise + 
                   Slepgrp + Spi + Fru + Cofe + Educgrp + Gretea + Menopause, 
                 data = MData_fem)

tidy(SurvM2, exponentiate = TRUE, conf.int = TRUE) %>% 
  knitr::kable(.)

