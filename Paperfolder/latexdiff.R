##%######################################################%##
#                                                          #
####                  make traced pdf                   ####
#                                                          #
##%######################################################%##
# library(latexdiffr)
# setwd("/Paperfolder")
# latexdiff("Paperfolder_v0_0_0.Rmd", "Paperfolder.Rmd", 
#           output = "diff_v0_0_1", clean = FALSE)
system("latexdiff Paperfolder_v0_0_0.tex Paperfolder.tex > diff_v_0_0_1.tex")

# abstract diff: Paperfolder_v0_0_0.tex Paperfolder.tex 20200716
\DIFaddbegin \DIFadd{In conclusion, %DIF >
    data from the JACC study has provided strong evidence that daily milk %DIF >
    intake among Japanese men was associated with delayed and lower hazard %DIF >
    of mortality from stroke especially cerebral infarction. }\DIFaddend

# calculate mortality rate using poisson regression
MData_men$Tot_Stro_bin <- MData_men$Tot_Stroke == "I60_9"
meanAge <- mean(MData_men$Age)
medianFY <- mean(MData_men$followpy)
MData_fem$Tot_Stro_bin <- MData_fem$Tot_Stroke == "I60_9"
meanAge <- mean(MData_fem$Age)
medianFY <- mean(MData_fem$followpy)




library(tidyverse)
ex <- MData_men %>% 
  select(Tot_Stro_bin, Age, Mlkfre, followpy)
ex <- MData_fem %>% 
  select(Tot_Stro_bin, Age, Mlkfre, followpy)

M <- glm(Tot_Stro_bin ~ Age + Mlkfre + offset(log(followpy)), data = ex)
# M <- glm(Tot_Stro_bin ~ Age + Mlkfre, data = ex)
summary(M)

newdata <- data.frame(Age = rep(meanAge, 5), 
                      Mlkfre = c( "Never", "Mon1_2", 
                                  "Wek1_2", "Wek3_4", "Daily"), 
                      followpy = medianFY)
# newdata <- data.frame(Age = rep(meanAge, 5), 
#                       Mlkfre = c( "Never", "Mon1_2", 
#                                   "Wek1_2", "Wek3_4", "Daily"))
Mpred <- predict(M, newdata, type="response",  se.fit = TRUE)
Mpred <- predict(M, newdata, type="link")
# lci  <-  Mpred$fit + (0.96 * Mpred$se.fit) 
# hci  <-  Mpred$fit - (0.95 * Mpred$se.fit)


MData_men %>%
  group_by(Mlkfre) %>%
  summarise(pyear = sum(followpy), n = n(), 
            N_stroke = sum(Tot_Stro_bin), 
            Rate = sum(Tot_Stro_bin)/sum(followpy) * 1000)
# 0.1773524 0.2024918 0.1717266 0.1610058 0.1517254 
# 0.1254686 0.1420398 0.1245280 0.1090470 0.1216790 

