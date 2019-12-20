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
#'   - \usepackage[ipa]{zxjafont} 
#' ---
#' # Read in the data

library(readr)
library(tidyverse)
library(lubridate) # for dealing with date time data 

MILK <- read_csv("data/StrokeMilk.csv", 
                     progress = show_progress(), 
                     col_types = cols(.default = "c"))

MILK %>% 
  filter(tr_age > 39 & tr_age < 80) %>% 
  group_by(tr_sex) %>% 
  summarise(n= n()) %>%
  mutate(rel.freq = paste0(round(100 * n/sum(n), 2), "%"))  

#' # delete subjects outside of age range ------------------------------------ 

MILK <- MILK %>%
  filter(tr_age > 39 & tr_age < 80)

#' # define total stroke mortality --------------------------------


# MILK <- 
  
MILK %>% 
  mutate(Tot_Stroke = if_else(grepl("I6[0-9][0-9]|I6[0-9]",  
                                    ICD10), "I60_9", 
                        if_else(!is.na(ICD10), "other_death", 
                                      "Alive/Censor"))) %>% 
  group_by(tr_sex, Tot_Stroke) %>% 
  summarise(n= n()) %>%
  mutate(rel.freq = paste0(round(100 * n/sum(n), 2), "%"))  



