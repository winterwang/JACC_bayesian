# jaccbayes


# stage <- larynx.dat$stage
Mlkfre <- as.numeric(MData_men$Mlkfre)
table(Mlkfre)
table(MData_men$Mlkfre) #1 = Never ; 2 = Mon1_2 ; 3 = Wek1_2 ; 4 = Wek3_4 ; 5 = Daily
MlkLogi <- as.numeric(MData_men$MlkLogi == "Drinker")

# is.censored <- is.na(t)
is.censored <- MData_men$Tot_Stroke != "I60_9"

# t <- larynx.dat$time
t <- if_else(!is.censored, MData_men$followpy, 0)
t <- na_if(t, 0)
t

# age <- larynx.dat$age
Agegrp <- as.numeric(MData_men$Agegrp)
table(Agegrp)
table(MData_men$Agegrp)  # 1 = [30,45) ; 2 = [45,55) ; 3 = [55,65) ; 4 =  [65,75) ; 5 = [75,80) 
Age <- MData_men$Age


# c <- larynx.dat$cens_time
c <- if_else(is.censored, MData_men$followpy, 0)


c[!is.censored] <- t[!is.censored] + 1
# c
# t
# is.censored


larynx.weibull.model <- function() {
  for(i in 1:39386){
    # sAge[i] <- (Age[i] - mean(Age[])) / sd(Age[])
    # sYr[i] <- (yr[i] - mean(yr[])) / sd(yr[])
    is.censored[i] ~ dinterval(t[i], c[i])
    t[i] ~ dweib(shape, lambda[i])
    lambda[i] <- exp(-mu[i] * shape)
    mu[i] <- beta[1] + beta[2]*equals(Mlkfre[i], 2) +
      beta[3]*equals(Mlkfre[i], 3) + beta[4] * equals(Mlkfre[i], 4) +
      beta[5] * equals(Mlkfre[i], 5) 
    
    ### calculate log-likelihoods
    # y[i] <- ifelse(is.censored[i], c[i], t[i])
    # loglik[i] <- log(ifelse(is.censored[i],
    #                         exp(-lambda[i] * (y[i] ^ shape)),
    #                         shape * lambda[i] * (y[i] ^ (shape - 1)) * exp(-lambda[i] * (y[i] ^ shape))))
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

# larynx.data <- c("t", "age", "yr", "c", "stage", "is.censored")
MILKdataMEN <- c("t",  "c", "Mlkfre", "is.censored")
larynx.params <- c("beta", "AFT", "HR", "p.crit", "shape")
library(R2jags)

# JACC.weibull.fit <- jags(data = MILKdataMEN,
#                          parameters.to.save = larynx.params,
#                          n.chains = 2,
#                          n.iter = 1000,
#                          n.burnin = 100,
#                          n.thin = 1,
#                          model.file = larynx.weibull.model)

jagsfit <- jags.parallel(data=MILKdataMEN,  parameters.to.save = larynx.params, 
                         n.iter=10000, n.burnin=(5000/2), n.chains = 3,
                         model.file=larynx.weibull.model)

# print(JACC.weibull.fit)
print(jagsfit)
# mcmcplots::traplot(JACC.weibull.fit, c("beta[5]", "HR[5]"))
mcmcplots::traplot(jagsfit, c("beta[5]", "HR[5]"))

# samplesHistory("*", mfrow = c(3,1), beg = 501, ask = FALSE)
Simulated <- coda::as.mcmc(JACC.weibull.fit)
library(ggmcmc)

ggSample <- ggs(Simulated)
ggSample %>% 
  filter(Iteration >= 500 & Parameter %in% c("beta[5]", "HR[5]")) %>% 
  ggs_traceplot()

# update(JACC.weibull.fit, n.iter = 1000)
# autojags(JACC.weibull.fit, n.iter = 5000, Rhat = 1.05, n.update = 2)

# postsamples <- buildMCMC("*")
# gelman.diag(Simulated)
# gelman.plot(Simulated)
Simulated <- coda::as.mcmc(JACC.weibull.fit)

ggSample <- ggs(Simulated)
ggSample %>% 
  filter(Iteration >= 100 & Parameter %in% c("beta[5]", "AFT[5]")) %>% 
  ggs_autocorrelation()

HR5 <- ggSample %>% 
  filter(Parameter == "HR[5]")

plot(density(HR5$value), main = "HR sample 10000", 
     ylab = "P(theta)", xlab = "theta", col = "red")
# hist(Y$value, main = "y sample 10000", ylab = "P(Y)", 
#      xlab = "y", col = "red", prob = TRUE)

##%######################################################%##
#                                                          #
####             age adjusted model in MEN              ####
#                                                          #
##%######################################################%##

Age.weibull.model <- function() {
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
jagsfit <- jags.parallel(data=MILKdataMEN,  parameters.to.save = age.params, 
                         n.iter=10000, n.burnin=(5000/2), n.chains = 3,
                         model.file=Age.weibull.model)
print(jagsfit)
# mcmcplots::traplot(JACC.weibull.fit, c("beta[5]", "HR[5]"))
mcmcplots::traplot(jagsfit, c("AFT[5]", "HR[5]"))

##%######################################################%##
#                                                          #
####       Model 3, smking, alcohol, bmi, DM, HT,       ####
####  KID, LIV, Excercise, Slep, Spi, Cofe, Educ, Gret  ####
#                                                          #
##%######################################################%##


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


M2.weibull.model <- function() {
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
jagsfit <- jags.parallel(data=MILKdataMEN,  parameters.to.save = M2.params, 
                         n.iter=10000, n.burnin=(5000/2), n.chains = 3,
                         model.file=M2.weibull.model)
print(jagsfit)
