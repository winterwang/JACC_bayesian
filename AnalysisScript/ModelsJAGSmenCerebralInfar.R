# HemoStroke == "I63"

library(tidyverse)
load("data/JACCmilkstrokewithHemo.Rdata")
#######################Model 0 crude#######################
#                                                          #
#                      Model 0 crude                       #
#                                                          #
#**********************************************************#


# Define exposure
Mlkfre <- as.numeric(MData_men$Mlkfre)
table(Mlkfre)
table(MData_men$Mlkfre) #1 = Never ; 2 = Mon1_2 ; 3 = Wek1_2 ; 4 = Wek3_4 ; 5 = Daily

# define events 
is.censored <- MData_men$HemoStroke != "I63"

# define followup time for events
t <- if_else(!is.censored, MData_men$followpy, 0)
t <- na_if(t, 0)
t


# define followup time for censored 
c <- if_else(is.censored, MData_men$followpy, 0)

c[!is.censored] <- t[!is.censored] + 1




jaccIschemic.weibull.model0 <- function() {
  for(i in 1:39386){
    is.censored[i] ~ dinterval(t[i], c[i])
    t[i] ~ dweib(shape, lambda[i])
    lambda[i] <- exp(-mu[i] * shape)
    mu[i] <- beta[1] + beta[2]*equals(Mlkfre[i], 2) +
      beta[3]*equals(Mlkfre[i], 3) + beta[4] * equals(Mlkfre[i], 4) +
      beta[5] * equals(Mlkfre[i], 5) 
  }
  
  ## priors for betas
  for(j in 1:5){
    beta[j] ~ dnorm(0, 0.001)
  }
  
  ### prior for shape
  shape ~ dgamma(.001, .001)
  
  ### Generated values
  AFT[2] <- exp(-beta[2])
  HR[2] <- exp(-shape * beta[2])
  p.crit[2] <- step(1 - HR[2])
  
  AFT[3] <- exp(-beta[3])
  HR[3] <- exp(-shape * beta[3])
  p.crit[3] <- step(1 - HR[3])
  
  AFT[4] <- exp(-beta[4])
  HR[4] <- exp(-shape * beta[4])
  p.crit[4] <- step(1 - HR[4])
  
  AFT[5] <- exp(-beta[5])
  HR[5] <- exp(-shape * beta[5])
  p.crit[5] <- step(1 - HR[5])
}

MILKdataMEN <- c("t",  "c", "Mlkfre", "is.censored")
m0.params <- c("beta", "AFT", "HR", "p.crit", "shape")
library(R2jags)


start.time <- Sys.time()
jagsfit <- jags.parallel(data=MILKdataMEN,  parameters.to.save = m0.params, 
                         n.iter=100000, n.burnin=(5000/2), n.chains = 4, #n.thin = 100,
                         model.file=jaccIschemic.weibull.model0)
end.time <- Sys.time()
end.time - start.time #Time difference of 9.715389 hours
M0menIscheStroke_20200702 <- jagsfit
print(M0menIscheStroke_20200702)
summary(mcmcplots::as.mcmc.rjags(M0menIscheStroke_20200702))

save(M0menIscheStroke_20200702, file = "IscheStrokeMen.RData")




#######################Model 1 age-adj######################
#                                                          #
#                      Model 1 age-adj                     #
#                                                          #
#**********************************************************#

Age <- MData_men$Age


AgeIsche.weibull.model <- function() {
  for(i in 1:39386){
    sAge[i] <- (Age[i] - mean(Age[])) / sd(Age[])
    # sYr[i] <- (yr[i] - mean(yr[])) / sd(yr[])
    is.censored[i] ~ dinterval(t[i], c[i])
    t[i] ~ dweib(shape, lambda[i])
    lambda[i] <- exp(-mu[i] * shape)
    mu[i] <- beta[1] + beta[2]*equals(Mlkfre[i], 2) +
      beta[3]*equals(Mlkfre[i], 3) + beta[4] * equals(Mlkfre[i], 4) +
      beta[5] * equals(Mlkfre[i], 5) + beta[6] * sAge[i]
  }
  
  ## priors for betas
  for(j in 1:6){
    beta[j] ~ dnorm(0, 0.001)
  }
  
  ### prior for shape
  shape ~ dgamma(.001, .001)
  
  ### Generated values
  AFT[2] <- exp(-beta[2])
  HR[2] <- exp(-shape * beta[2])
  p.crit[2] <- step(1 - HR[2])
  
  AFT[3] <- exp(-beta[3])
  HR[3] <- exp(-shape * beta[3])
  p.crit[3] <- step(1 - HR[3])
  
  AFT[4] <- exp(-beta[4])
  HR[4] <- exp(-shape * beta[4])
  p.crit[4] <- step(1 - HR[4])
  
  AFT[5] <- exp(-beta[5])
  HR[5] <- exp(-shape * beta[5])
  p.crit[5] <- step(1 - HR[5])
  
}

MILKdataMEN <- c("t",  "c", "Mlkfre", "is.censored", "Age")
age.params <- c("AFT", "HR", "p.crit")
start.time <- Sys.time()
jagsfit <- jags.parallel(data=MILKdataMEN,  parameters.to.save = age.params, 
                         n.iter=100000, n.burnin=(5000/2), n.chains = 4,
                         model.file=AgeIsche.weibull.model)
end.time <- Sys.time()
end.time - start.time #

M1menIscheStroke_20200702 <- jagsfit
print(M1menIscheStroke_20200702)
summary(mcmcplots::as.mcmc.rjags(M1menIscheStroke_20200702))
save(M0menIscheStroke_20200702, M1menIscheStroke_20200702, file = "IscheStrokeMen.RData")




#######################Model 2 Mult-ad######################
#                                                          #
#                      Model 2 Mult-ad                     #
#                                                          #
#**********************************************************#



Smoking <- as.numeric(MData_men$Smoking)
table(Smoking)
table(MData_men$Smoking) #1 = Never ; 2 = Past ; 3 = Current ; 4 = unknown 

Alc_Fre <- as.numeric(MData_men$Alc_Fre)
table(Alc_Fre)
table(MData_men$Alc_Fre) #1 = < 1/week ; 2 = 1-4 /week ; 3 = Daily ; 4 = Never or past; 5 = unknown

BMI <- as.numeric(MData_men$BMI)
BMIgrp <- as.numeric(MData_men$BMIgrp)
table(BMIgrp)
table(MData_men$BMIgrp) #1 = [18.5,25) ; 2 = [14,18.5) ; 3 = [25,30) ; 4 = [30,40); 5 = unknown

DM_hist <- as.numeric(as.factor(MData_men$DM_hist))
table(DM_hist)
table(MData_men$DM_hist) #1 =  FALSE ; 2 = TRUE ; 3 = unknown

HT_hist <- as.numeric(as.factor(MData_men$HT_hist))
table(HT_hist)
table(MData_men$HT_hist) #1 =  FALSE ; 2 = TRUE ; 3 = unknown

KID_hist <- as.numeric(as.factor(MData_men$KID_hist))
table(KID_hist)
table(MData_men$KID_hist) #1 =  FALSE ; 2 = TRUE ; 3 = unknown


LIV_hist <- as.numeric(as.factor(MData_men$LIV_hist))
table(LIV_hist)
table(MData_men$LIV_hist) #1 =  FALSE ; 2 = TRUE ; 3 = unknown

Exercise <- as.numeric(MData_men$Exercise)
table(Exercise)
table(MData_men$Exercise) #1 =  Almost0  ; 2 = > 1h/w ; 3 = unknown

Slepgrp <- as.numeric(MData_men$Slepgrp)
table(Slepgrp)
table(MData_men$Slepgrp) #1 =  [0,6.9)  ; 2 = [6.9,7.9) ; 3 = [7.9,8.9) ; 4 = [8.9,23); 5 = unknown

Spi <- as.numeric(MData_men$Spi)
table(Spi)
table(MData_men$Spi) #1 =  Less1tm  ; 2 = One2tw ; 3 = Thre4tw ; 4 = daily ; 5 = unknown

Fru <- as.numeric(MData_men$Fru)
table(Fru)
table(MData_men$Fru) #1 =  Less1tm  ; 2 = One2tw ; 3 = Thre4tw ; 4 = daily ; 5 = unknown


Cofe <- as.numeric(MData_men$Cofe)
table(Cofe)
table(MData_men$Cofe) #1 =  daily  ; 2 = Thre3tw  ; 3 = Never ; 4 = unknown


Educgrp <- as.numeric(MData_men$Educgrp)
table(Educgrp)
table(MData_men$Educgrp) #1 =  [0,18)  ; 2 = [18,70)  ; 3 = unknown

Gretea <- as.numeric(MData_men$Gretea)
table(Gretea)
table(MData_men$Gretea) #1 =  daily   ; 2 = Thre3tw  ; 3 = Never ; 4 = unknown




M2hemo.weibull.model <- function() {
  for(i in 1:39386){
    sAge[i] <- (Age[i] - mean(Age[])) / sd(Age[])
    # sBMI[i] <- (BMI[i] - mean(BMI[])) / sd(BMI[])
    is.censored[i] ~ dinterval(t[i], c[i])
    t[i] ~ dweib(shape, lambda[i])
    lambda[i] <- exp(-mu[i] * shape)
    mu[i] <- beta[1] + beta[2]*equals(Mlkfre[i], 2) +
      beta[3]*equals(Mlkfre[i], 3) + beta[4] * equals(Mlkfre[i], 4) +
      beta[5] * equals(Mlkfre[i], 5) + beta[6] * sAge[i] + beta[7] * equals(Smoking[i], 2) + 
      beta[8] * equals(Smoking[i], 3) + beta[9] * equals(Smoking[i], 4) + 
      beta[10] * equals(Alc_Fre[i], 1) + beta[11] * equals(Alc_Fre[i], 2) + 
      beta[12] * equals(Alc_Fre[i], 3) + beta[13] * equals(Alc_Fre[i], 5) + 
      beta[14] * equals(BMIgrp[i], 2) + beta[15] * equals(BMIgrp[i], 3) + 
      beta[16] * equals(BMIgrp[i], 4) + beta[17] * equals(BMIgrp[i], 5)
    #+ beta[15] * equals(DM_hist[i], 2) +  
    # beta[16] * equals(DM_hist[i], 3) + 
    # beta[17] * equals(HT_hist[i], 2) + beta[18] * equals(HT_hist[i], 3) + 
    # beta[19] * equals(KID_hist[i], 2) + beta[20] * equals(KID_hist[i], 3)  + 
    # beta[21] * equals(LIV_hist[i], 2) + beta[22] * equals(LIV_hist[i], 3)  + 
    # beta[23] * equals(Exercise[i], 2) + beta[24] * equals(Exercise[i], 3)  +
    # beta[25] * equals(Slepgrp[i], 2) + beta[26] * equals(Slepgrp[i], 3) + 
    # beta[27] * equals(Slepgrp[i], 4) + beta[28] * equals(Slepgrp[i], 5) +
    # beta[29] * equals(Spi[i], 2) + beta[30] * equals(Spi[i], 3) + 
    # beta[31] * equals(Spi[i], 4) + beta[32] * equals(Spi[i], 5) + 
    # beta[33] * equals(Cofe[i], 2) + beta[34] * equals(Cofe[i], 3) + 
    # beta[35] * equals(Cofe[i], 4) + beta[36] * equals(Educgrp[i], 2) + 
    # beta[37] * equals(Educgrp[i], 3) +  beta[38] * equals(Gretea[i], 2) + 
    # beta[39] * equals(Gretea[i], 3) +  beta[40] * equals(Gretea[i], 4) + 
    # beta[41] * equals(Fru[i], 2) + beta[42] * equals(Fru[i], 3) + 
    # beta[43] * equals(Fru[i], 4) + beta[44] * equals(Fru[i], 5) #+ 
    # beta[45] * equals(BMIgrp[i], 3) + beta[46] * equals(BMIgrp[i], 4) + 
    # beta[47] * equals(BMIgrp[i], 5) 
  }
  
  ## priors for betas
  for(j in 1:17){
    beta[j] ~ dnorm(0, 0.001)
  }
  
  ### prior for shape
  shape ~ dgamma(.001, .001)
  
  ### Generated values
  AFT[2] <- exp(-beta[2])
  HR[2] <- exp(-shape * beta[2])
  p.crit[2] <- step(1 - HR[2])
  
  AFT[3] <- exp(-beta[3])
  HR[3] <- exp(-shape * beta[3])
  p.crit[3] <- step(1 - HR[3])
  
  AFT[4] <- exp(-beta[4])
  HR[4] <- exp(-shape * beta[4])
  p.crit[4] <- step(1 - HR[4])
  
  AFT[5] <- exp(-beta[5])
  HR[5] <- exp(-shape * beta[5])
  p.crit[5] <- step(1 - HR[5])
  
}

MILKdataMEN <- c("t",  "c", "Mlkfre", "is.censored", "Age", "Smoking", 
                 "Alc_Fre", "BMIgrp"#, 
                 # "DM_hist", "HT_hist", "KID_hist", 
                 # "LIV_hist", "Exercise", "Slepgrp", "Spi", "Cofe", 
                 # "Educgrp", "Gretea", "Fru"
)
M2.params <- c("AFT", "HR", "p.crit")
start.time <- Sys.time()
jagsfit <- jags.parallel(data=MILKdataMEN,  parameters.to.save = M2.params, 
                         n.iter=100000, n.burnin=(3000/2), n.chains = 4,
                         model.file=M2hemo.weibull.model)
end.time <- Sys.time()
end.time - start.time 

M2menIscheStroke_20200702 <- jagsfit
print(M2menIscheStroke_20200702)
summary(mcmcplots::as.mcmc.rjags(M2menIscheStroke_20200702))
save(M0menIscheStroke_20200702, 
     M1menIscheStroke_20200702, 
     M2menIscheStroke_20200702, file = "IscheStrokeMen.RData")


