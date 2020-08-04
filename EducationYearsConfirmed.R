MILK <- read_csv("data/StrokeMilk.csv", 
                 progress = show_progress(), 
                 col_types = cols(.default = "c"))


MILK_0 <- MILK %>%
  filter(tr_age > 39 & tr_age < 80)



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




a <- MILK_0 %>% 
  mutate(Educ = as.numeric(MILK_0$SCHOOL)) %>% 
  mutate(Educgrp = cut(Educ, breaks = c(0, 19, 70), right = FALSE)) %>% 
  mutate(Educgrp = as.character(Educgrp)) %>% 
  replace_na(list(Educgrp = "unknown")) 



a %>% 
  group_by(tr_sex,  Milk_fre, Educgrp) %>% 
  summarise (n= n()) %>%
  mutate(rel.freq = paste0(round(100 * n/sum(n), 2), "%"))  %>% 
  print(n=Inf)

a %>% 
  group_by(tr_sex,  MlkLogi, Educgrp) %>% 
  summarise (n= n()) %>%
  mutate(rel.freq = paste0(round(100 * n/sum(n), 2), "%"))  %>% 
  print(n=Inf)


a %>% 
  mutate(ENEy = as.numeric(ENERGY)) %>% 
  group_by(tr_sex,  Milk_fre) %>% 
  summarise (MeanEnery = mean(ENEy, na.rm = T), SDEnergy = sd(ENEy, na.rm = T)) %>%
  # mutate(rel.freq = paste0(round(100 * n/sum(n), 2), "%"))  %>% 
  print(n=Inf)

a %>% 
  mutate(ENEy = as.numeric(ENERGY)) %>% 
  group_by(tr_sex,  MlkLogi) %>% 
  summarise (MeanEnery = mean(ENEy, na.rm = T), SDEnergy = sd(ENEy, na.rm = T)) %>%
  # mutate(rel.freq = paste0(round(100 * n/sum(n), 2), "%"))  %>% 
  print(n=Inf)

