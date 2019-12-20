##%######################################################%##
#                                                          #
####                    data how to                     ####
#                                                          #
##%######################################################%##

##1. 127208 subjects at first 53931 male and 73277 female 
##2. 110585 aged between 40-79: 46395 male and 64190 female 
Milk <- subset(Milk_dataset, (Age < 80)&(Age >= 40))
##3. Previous history of cancer or stroke or chronic CVD 
#####(myocardial infarction, angina pectoris, and other chonic 
#####ischemic heart diseases) n = 5693 ====>> excluded
use(Milk)
length(ID[(Canhistory1 == "Yes")|(Canhistory2 == "Yes")]) #n = 1461
length(ID[(Apohistory == "Yes")]) #n = 1496
length(ID[(Canhistory1 == "Yes")|(Canhistory2 == "Yes")|
            (Apohistory == "Yes")|(MIhistory == "Yes")|
            (Anginahistory == "Yes")|(OtherCirhistory == 1)]) 
Milk_ill_excluded <- subset(Milk, (Canhistory1 != "Yes")&(Canhistory2 != "Yes")&
                              (Apohistory != "Yes")&(MIhistory != "Yes")&
                              (Anginahistory != "Yes")&(OtherCirhistory != 1))
##### 104892 subjects 43902 male and 60990 female 

##4. data missing on milk intake frequency n = 9912, 4263 men and 5649 women ==> excluded 
MILK_ill_milk_excluded <- subset(Milk_ill_excluded, 
                                 Milk_ill_excluded$Milk != "unknown")

##FINALLY 
tab1(MILK_ill_milk_excluded$Gender)
### 94980 subjects (39639 male and 55341 female)
use(MILK_ill_milk_excluded)
MILK_ill_milk_excluded$Milk<- MILK_ill_milk_excluded$Milk[drop = TRUE]




##4. Death within 5 years of follow-up n = 4922 ==> excluded 
use(MILK_ill_milk_excluded)
summ(Follow_year)
length(ID[Follow_year <= 5]) ##n = 4922
  Milk_illandearlydeath_excluded <- subset(MILK_ill_milk_excluded, Follow_year > 5)
  ##### 90058, 37002 men and 53056 women
tab1(MILK_ill_milk_excluded$Milk)


##%######################################################%##
#                                                          #
####                   data cleaning                    ####
#                                                          #
##%######################################################%##


# etwd("~/Dropbox/JACC _study/dataset")
library(epicalc)
Call01[Call01==". "] <- NA
Call02[Call02==". "] <- NA
Call03[Call03==". "] <- NA
Call03[Call03=="X "] <- NA
Call04[Call04==". "] <- NA
Call05[Call05==". "] <- NA
Call05[Call05=="X "] <- NA
Call06[Call06 == '. '] <- NA
Call06[Call06 == 'X '] <- NA

Areano_region <- as.integer(Call01$areano)
Gender <- as.factor(Call01$sex)
Age <- as.numeric(Call01$age)
School <- as.numeric(Call06$school)

height <- as.numeric(Call06$ht10)
weight <- as.numeric(Call06$wt10)
Body.Mass.Index <- (weight/10)/((height/1000)^2)
Sport <- as.factor(Call03$SPORT)
tab1(Sport)

Sampo <- as.factor(Call03$SAMPO)
tab1(Sampo)



### disease histories 
Apohistory <- as.character(Call01$p_apo)
Apohistory[as.numeric(Apohistory) > 1] <- "Yes"
Apohistory[is.na(Apohistory)] <- "unknown"
Apohistory[as.numeric(Apohistory) == 1] <- "No"
Apohistory <- as.factor(Apohistory)
print(levels(Apohistory))
Apohistory <- factor(Apohistory, levels(Apohistory)[c(3,1,2)])

Hythistory <- as.character(Call01$p_ht)
tab1(Hythistory)
Hythistory[as.numeric(Hythistory) > 1] <- "Yes"
Hythistory[is.na(Hythistory)] <- "unknown"
Hythistory[as.numeric(Hythistory) == 1] <- "No"
Hythistory <- as.factor(Hythistory)
print(levels(Hythistory))
Hythistory <- factor(Hythistory, levels(Hythistory)[c(3,1,2)])


MIhistory <- as.character(Call01$p_mi)
tab1(MIhistory)
MIhistory[as.numeric(MIhistory) > 1] <- "Yes"
MIhistory[is.na(MIhistory)] <- "unknown"
MIhistory[as.numeric(MIhistory) == 1] <- "No"
MIhistory <- as.factor(MIhistory)
print(levels(MIhistory))
MIhistory <- factor(MIhistory, levels(MIhistory)[c(3,1,2)])

Livhistory <- as.character(Call01$p_liv)
Livhistory[Livhistory == "X "] <- NA
tab1(Livhistory)
Livhistory[as.numeric(Livhistory) > 1] <- "Yes"
Livhistory[is.na(Livhistory)] <- "unknown"
Livhistory[as.numeric(Livhistory) == 1] <- "No"
Livhistory <- as.factor(Livhistory)
print(levels(Livhistory))
Livhistory <- factor(Livhistory, levels(Livhistory)[c(3,1,2)])


DMhistory <- as.character(Call01$p_dm)
tab1(DMhistory)
DMhistory[as.numeric(DMhistory) > 1] <- "Yes"
DMhistory[is.na(DMhistory)] <- "unknown"
DMhistory[as.numeric(DMhistory) == 1] <- "No"
DMhistory <- as.factor(DMhistory)
print(levels(DMhistory))
DMhistory <- factor(DMhistory, levels(DMhistory)[c(3,1,2)])


Canhistory1 <- as.character(Call01$p_can1)
Canhistory1[Canhistory1 == "X "] <- NA
tab1(Canhistory1)
Canhistory1[as.numeric(Canhistory1) > 1] <- "Yes"
Canhistory1[is.na(Canhistory1)] <- "unknown"
Canhistory1[as.numeric(Canhistory1) == 1] <- "No"
Canhistory1 <- as.factor(Canhistory1)
print(levels(Canhistory1))
Canhistory1 <- factor(Canhistory1, levels(Canhistory1)[c(3,1,2)])

Canhistory2 <- as.character(Call01$p_can2)
Canhistory2[Canhistory2 == "X "] <- NA
tab1(Canhistory2)
Canhistory2[as.numeric(Canhistory2) > 1] <- "Yes"
Canhistory2[is.na(Canhistory2)] <- "unknown"
Canhistory2[as.numeric(Canhistory2) == 1] <- "No"
Canhistory2 <- as.factor(Canhistory2)
print(levels(Canhistory2))
Canhistory2 <- factor(Canhistory2, levels(Canhistory2)[c(3,1,2)])

Anginahistory <- 1:127208
Anginahistory[(as.numeric(Call01$p_oth1x) == 413)|
                (as.numeric(Call01$p_oth2x) == 413)] <- TRUE
Anginahistory[Anginahistory != TRUE] <- FALSE
tab1(Anginahistory)

OtherCirhistory <- 1:127208
OtherCirhistory[(as.numeric(Call01$p_oth1x) == 414)|
                  (as.numeric(Call01$p_oth2x) == 414)] <- TRUE
OtherCirhistory[OtherCirhistory != TRUE] <- FALSE
tab1(OtherCirhistory)


ID  <- 1:127208
dataset <- data.frame(ID, Areano_region, Age, Gender, weight, 
                      height, Body.Mass.Index, School,
                      Sport, Sampo, Anginahistory,
                      Apohistory, Canhistory1, Canhistory2, DMhistory,
                      Hythistory, Livhistory, MIhistory, OtherCirhistory)

use(dataset)
tab1(Gender)
#Gender :  
#       Frequency Percent Cum. percent
#1men       53931    42.4         42.4
#2women     73277    57.6        100.0
#Total     127208   100.0        100.0
JACC_confirm <- subset(dataset, (Age > 39)&(Age < 80))
#n = 110585 
use(JACC_confirm)
tab1(JACC_confirm$Gender)
#Gender :  
#        Frequency Percent Cum. percent
#1men        46395    42            42
#2women      64190    58         100.0
#Total      110585   100.0       100.0

JACC_confirm0 <- subset(JACC_confirm, )

tab1((JACC_confirm$Apohistory == "Yes")|(JACC_confirm$MIhistory == "Yes"))
tab1((JACC_confirm$Apohistory == "Yes")|(JACC_confirm$MIhistory == "Yes")|
       (JACC_confirm$Anginahistory == 1)|(JACC_confirm$OtherCirhistory == 1))
tab1(JACC_confirm$OtherCirhistory)


tab1(JACC_confirm$School)

tab1(JACC_confirm$Sport)
tab1(JACC_confirm$Sampo)
tab1(is.na(JACC_confirm$Sport)|(is.na(JACC_confirm$Sampo)))

summary(JACC_confirm$Body.Mass.Index)


##%######################################################%##
#                                                          #
####               mortality from sroke                 ####
#                                                          #
##%######################################################%##

#### find people died from total stroke "ICD10 code: I60-I69"
use(MILK_ill_milk_excluded)
Stroke <- grep("^I6",MILK_ill_milk_excluded$ICD10)
MILK_ill_milk_excluded$siin <- rep("NO",nrow(MILK_ill_milk_excluded))
MILK_ill_milk_excluded$siin[Stroke] <- "Stroke"
table(MILK_ill_milk_excluded$siin)
MILK_ill_milk_excluded$DeathStroke[siin == "Stroke"] <- TRUE
table(MILK_ill_milk_excluded$DeathStroke)


