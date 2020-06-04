##%######################################################%##
#                                                          #
####                  data preparation                  ####
#                                                          #
##%######################################################%##

# stage <- larynx.dat$stage
Mlkfre <- as.numeric(MData_fem$Mlkfre)
table(Mlkfre)
table(MData_fem$Mlkfre) #1 = Never ; 2 = Mon1_2 ; 3 = Wek1_2 ; 4 = Wek3_4 ; 5 = Daily
MlkLogi <- as.numeric(MData_fem$MlkLogi == "Drinker")

# is.censored <- is.na(t)
is.censored <- MData_fem$Tot_Stroke != "I60_9"

# t <- larynx.dat$time
t <- if_else(!is.censored, MData_fem$followpy, 0)
t <- na_if(t, 0)
t

# age <- larynx.dat$age
Agegrp <- as.numeric(MData_fem$Agegrp)
table(Agegrp)
table(MData_fem$Agegrp)  # 1 = [30,45) ; 2 = [45,55) ; 3 = [55,65) ; 4 =  [65,75) ; 5 = [75,80) 
Age <- MData_fem$Age


# c <- larynx.dat$cens_time
c <- if_else(is.censored, MData_fem$followpy, 0)


c[!is.censored] <- t[!is.censored] + 1
# c
# t
# is.censored


larynx.weibull.model <- function() {
  for(i in 1:54999){
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

MILKdatafem <- c("t",  "c", "Mlkfre", "is.censored")
larynx.params <- c("beta", "AFT", "HR", "p.crit", "shape")
library(R2jags)

start.time <- Sys.time()
jagsfit <- jags.parallel(data=MILKdatafem,  parameters.to.save = larynx.params, 
                         n.iter=100000, n.burnin=(5000/2), n.chains = 3,
                         model.file=larynx.weibull.model)
end.time <- Sys.time()
end.time - start.time #Time difference of 1.316827  days
M0fem_20200526 <- jagsfit
# print(JACC.weibull.fit)
print(M0fem_20200526)
# Inference for Bugs model at "larynx.weibull.model", fit using jags,
# 3 chains, each with 1e+05 iterations (first 2500 discarded), n.thin = 97
# n.sims = 3015 iterations saved
#             mu.vect sd.vect      2.5%       25%       50%       75%     97.5%  Rhat n.eff
# AFT[2]        0.884   0.070     0.754     0.836     0.881     0.929     1.029 1.001  3000
# AFT[3]        0.874   0.052     0.778     0.840     0.872     0.908     0.978 1.002  1300
# AFT[4]        0.799   0.049     0.705     0.765     0.797     0.831     0.902 1.003  2200
# AFT[5]        0.877   0.041     0.800     0.850     0.876     0.903     0.958 1.007   470
# HR[2]         0.825   0.103     0.643     0.755     0.820     0.890     1.046 1.001  2900
# HR[3]         0.811   0.075     0.675     0.760     0.807     0.858     0.966 1.002  1100
# HR[4]         0.704   0.068     0.581     0.658     0.700     0.746     0.851 1.003  1700
# HR[5]         0.814   0.059     0.706     0.775     0.812     0.852     0.934 1.007   420
# beta[1]       5.107   0.074     4.985     5.065     5.104     5.144     5.242 1.030   210
# beta[2]       0.127   0.079    -0.028     0.074     0.126     0.179     0.283 1.001  3000
# beta[3]       0.136   0.059     0.022     0.097     0.137     0.174     0.251 1.002  1300
# beta[4]       0.226   0.062     0.104     0.185     0.226     0.268     0.350 1.003  2200
# beta[5]       0.133   0.046     0.043     0.102     0.132     0.163     0.223 1.007   470
# p.crit[2]     0.946   0.227     0.000     1.000     1.000     1.000     1.000 1.005  1800
# p.crit[3]     0.987   0.114     1.000     1.000     1.000     1.000     1.000 1.034  1100
# p.crit[4]     0.999   0.026     1.000     1.000     1.000     1.000     1.000 1.291  1500
# p.crit[5]     0.996   0.066     1.000     1.000     1.000     1.000     1.000 1.215   380
# shape         1.572   0.041     1.491     1.548     1.573     1.597     1.647 1.010   320
# deviance  17083.697  72.790 16944.439 17034.564 17081.900 17131.130 17229.728 1.001  3000
# 
# For each parameter, n.eff is a crude measure of effective sample size,
# and Rhat is the potential scale reduction factor (at convergence, Rhat=1).
# 
# DIC info (using the rule, pD = var(deviance)/2)
# pD = 2649.7 and DIC = 19733.4
# DIC is an estimate of expected predictive error (lower deviance is better).

save.image(file = "data/JACCmilkstroke.Rdata")

##%######################################################%##
#                                                          #
####             age adjusted model in Women              ####
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

MILKdatafem <- c("t",  "c", "Mlkfre", "is.censored", "Age")
age.params <- c("AFT", "HR", "p.crit")
start.time <- Sys.time()
jagsfit <- jags.parallel(data=MILKdatafem,  parameters.to.save = age.params, 
                         n.iter=100000, n.burnin=(500/2), n.chains = 3,
                         model.file=Age.weibull.model)
end.time <- Sys.time() 
end.time - start.time  #Time difference of 15.50321 hours
#print(jagsfit)
M1fem_20200603 <- jagsfit
# print(JACC.weibull.fit)
print(M1fem_20200603)
# Inference for Bugs model at "Age.weibull.model", fit using jags,
# 3 chains, each with 1e+05 iterations (first 250 discarded), n.thin = 99
# n.sims = 3021 iterations saved
#             mu.vect sd.vect      2.5%       25%       50%       75%     97.5%  Rhat n.eff
# AFT[2]        0.999   0.088     0.847     0.942     0.994     1.051     1.172 1.003  1700
# AFT[3]        1.106   0.085     0.973     1.057     1.100     1.148     1.257 1.004   640
# AFT[4]        1.017   0.078     0.891     0.970     1.014     1.058     1.164 1.003  1200
# AFT[5]        0.950   0.058     0.862     0.917     0.946     0.978     1.055 1.004  1600
# HR[2]         1.001   0.143     0.757     0.903     0.990     1.088     1.309 1.002  1500
# HR[3]         1.189   0.138     0.953     1.100     1.179     1.268     1.479 1.004   610
# HR[4]         1.030   0.123     0.820     0.950     1.024     1.102     1.281 1.003   900
# HR[5]         0.916   0.086     0.775     0.861     0.908     0.962     1.092 1.004  1100
# p.crit[2]     0.523   0.500     0.000     0.000     1.000     1.000     1.000 1.002  1200
# p.crit[3]     0.063   0.243     0.000     0.000     0.000     0.000     1.000 1.006  1300
# p.crit[4]     0.420   0.494     0.000     0.000     0.000     1.000     1.000 1.002  1400
# p.crit[5]     0.868   0.339     0.000     1.000     1.000     1.000     1.000 1.007   580
# deviance  11991.178  79.834 11866.165 11945.083 11985.542 12028.868 12125.671 1.003  1300
# 
# For each parameter, n.eff is a crude measure of effective sample size,
# and Rhat is the potential scale reduction factor (at convergence, Rhat=1).
# 
# DIC info (using the rule, pD = var(deviance)/2)
# pD = 3184.1 and DIC = 15175.3
# DIC is an estimate of expected predictive error (lower deviance is better).
save.image(file = "data/JACCmilkstroke.Rdata")


# mcmcplots::traplot(JACC.weibull.fit, c("beta[5]", "HR[5]"))
mcmcplots::traplot(jagsfit, c("AFT[5]", "HR[5]"))




##%######################################################%##
#                                                          #
####       Model 3, smking, alcohol, bmi, DM, HT,       ####
####  KID, LIV, Excercise, Slep, Spi, Cofe, Educ, Gret  ####
#                                                          #
##%######################################################%##


Smoking <- as.numeric(MData_fem$Smoking)
table(Smoking)
table(MData_fem$Smoking) #1 = Never ; 2 = Past ; 3 = Current ; 4 = unknown 

Alc_Fre <- as.numeric(MData_fem$Alc_Fre)
table(Alc_Fre)
table(MData_fem$Alc_Fre) #1 = < 1/week ; 2 = 1-4 /week ; 3 = Daily ; 4 = Never or past; 5 = unknown

BMI <- as.numeric(MData_fem$BMI)
BMIgrp <- as.numeric(MData_fem$BMIgrp)
table(BMIgrp)
table(MData_fem$BMIgrp) #1 = [18.5,25) ; 2 = [14,18.5) ; 3 = [25,30) ; 4 = [30,40); 5 = unknown

DM_hist <- as.numeric(as.factor(MData_fem$DM_hist))
table(DM_hist)
table(MData_fem$DM_hist) #1 =  FALSE ; 2 = TRUE ; 3 = unknown

HT_hist <- as.numeric(as.factor(MData_fem$HT_hist))
table(HT_hist)
table(MData_fem$HT_hist) #1 =  FALSE ; 2 = TRUE ; 3 = unknown

KID_hist <- as.numeric(as.factor(MData_fem$KID_hist))
table(KID_hist)
table(MData_fem$KID_hist) #1 =  FALSE ; 2 = TRUE ; 3 = unknown


LIV_hist <- as.numeric(as.factor(MData_fem$LIV_hist))
table(LIV_hist)
table(MData_fem$LIV_hist) #1 =  FALSE ; 2 = TRUE ; 3 = unknown

Exercise <- as.numeric(MData_fem$Exercise)
table(Exercise)
table(MData_fem$Exercise) #1 =  Almost0  ; 2 = > 1h/w ; 3 = unknown

Slepgrp <- as.numeric(MData_fem$Slepgrp)
table(Slepgrp)
table(MData_fem$Slepgrp) #1 =  [0,6.9)  ; 2 = [6.9,7.9) ; 3 = [7.9,8.9) ; 4 = [8.9,23); 5 = unknown

Spi <- as.numeric(MData_fem$Spi)
table(Spi)
table(MData_fem$Spi) #1 =  Less1tm  ; 2 = One2tw ; 3 = Thre4tw ; 4 = daily ; 5 = unknown

Fru <- as.numeric(MData_fem$Fru)
table(Fru)
table(MData_fem$Fru) #1 =  Less1tm  ; 2 = One2tw ; 3 = Thre4tw ; 4 = daily ; 5 = unknown


Cofe <- as.numeric(MData_fem$Cofe)
table(Cofe)
table(MData_fem$Cofe) #1 =  daily  ; 2 = Thre3tw  ; 3 = Never ; 4 = unknown


Educgrp <- as.numeric(MData_fem$Educgrp)
table(Educgrp)
table(MData_fem$Educgrp) #1 =  [0,18)  ; 2 = [18,70)  ; 3 = unknown

Gretea <- as.numeric(MData_fem$Gretea)
table(Gretea)
table(MData_fem$Gretea) #1 =  daily   ; 2 = Thre3tw  ; 3 = Never ; 4 = unknown


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
      beta[16] * equals(BMIgrp[i], 4) + beta[17] * equals(BMIgrp[i], 5) +
      beta[18] * equals(DM_hist[i], 2) + beta[19] * equals(DM_hist[i], 3) + 
    beta[20] * equals(HT_hist[i], 2) + beta[21] * equals(HT_hist[i], 3) +
    beta[22] * equals(KID_hist[i], 2) + beta[23] * equals(KID_hist[i], 3)  +
    beta[24] * equals(LIV_hist[i], 2) + beta[25] * equals(LIV_hist[i], 3)  +
    beta[26] * equals(Exercise[i], 2) + beta[27] * equals(Exercise[i], 3)  +
    beta[28] * equals(Slepgrp[i], 2) + beta[29] * equals(Slepgrp[i], 3) +
    beta[30] * equals(Slepgrp[i], 4) + beta[31] * equals(Slepgrp[i], 5) +
    beta[32] * equals(Spi[i], 2) + beta[33] * equals(Spi[i], 3) +
    beta[34] * equals(Spi[i], 4) + beta[35] * equals(Spi[i], 5) +
    beta[36] * equals(Cofe[i], 2) + beta[37] * equals(Cofe[i], 3) +
    beta[38] * equals(Cofe[i], 4) + beta[39] * equals(Educgrp[i], 2) +
    beta[40] * equals(Educgrp[i], 3) +  beta[41] * equals(Gretea[i], 2) +
    beta[42] * equals(Gretea[i], 3) +  beta[43] * equals(Gretea[i], 4) +
    beta[44] * equals(Fru[i], 2) + beta[45] * equals(Fru[i], 3) +
    beta[46] * equals(Fru[i], 4) + beta[47] * equals(Fru[i], 5) #+
   
  }
  
  ## priors for betas
  for(j in 1:47){
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

MILKdatafem <- c("t",  "c", "Mlkfre", "is.censored", "Age", "Smoking", 
                 "Alc_Fre", "BMIgrp", 
                 "DM_hist", "HT_hist", "KID_hist",
                 "LIV_hist", "Exercise", "Slepgrp", "Spi", "Cofe",
                 "Educgrp", "Gretea", "Fru"
)
M2.params <- c("AFT", "HR", "p.crit")
jagsfit <- jags.parallel(data=MILKdatafem,  parameters.to.save = M2.params, 
                         n.iter=10000, n.burnin=(5000/2), n.chains = 3,
                         model.file=M2.weibull.model)
print(jagsfit)
