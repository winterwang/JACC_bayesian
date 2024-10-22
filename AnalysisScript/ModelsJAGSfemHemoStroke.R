# HemoStroke == "I60_2"

library(tidyverse)
load("data/JACCmilkstrokewithHemo.Rdata")
#######################Model 0 crude#######################
#                                                          #
#                      Model 0 crude                       #
#                                                          #
#**********************************************************#


# Define exposure
Mlkfre <- as.numeric(MData_fem$Mlkfre)
table(Mlkfre)
table(MData_fem$Mlkfre) #1 = Never ; 2 = Mon1_2 ; 3 = Wek1_2 ; 4 = Wek3_4 ; 5 = Daily

# define events 
is.censored <- MData_fem$HemoStroke != "I60_2"

# define followup time for events
t <- if_else(!is.censored, MData_fem$followpy, 0)
t <- na_if(t, 0)
t


# define followup time for censored 
c <- if_else(is.censored, MData_fem$followpy, 0)

c[!is.censored] <- t[!is.censored] + 1




jacchemo.weibull.model0 <- function() {
  for(i in 1:54999){
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

MILKdataFEM <- c("t",  "c", "Mlkfre", "is.censored")
m0.params <- c("beta", "AFT", "HR", "p.crit", "shape")
library(R2jags)


start.time <- Sys.time()
jagsfit <- jags.parallel(data=MILKdataFEM,  parameters.to.save = m0.params, 
                         n.iter=100000, n.burnin=(5000/2), n.chains = 4, #n.thin = 100,
                         model.file=jacchemo.weibull.model0)
end.time <- Sys.time()
end.time - start.time #Time difference of 22.07081 hours
M0FEMHemoStroke_20200704 <- jagsfit
print(M0FEMHemoStroke_20200704)
summary(mcmcplots::as.mcmc.rjags(M0FEMHemoStroke_20200704))
# save.image(file = "data/JACCmilkstrokewithHemo.Rdata")
# save.image(file = "data/JACCmilkstroke.Rdata")
save(M0FEMHemoStroke_20200704, file = "HemoStrokeFEM.RData")

# Inference for Bugs model at "jacchemo.weibull.model0", fit using jags,
# 4 chains, each with 1e+05 iterations (first 2500 discarded), n.thin = 97
# n.sims = 4020 iterations saved
#            mu.vect sd.vect     2.5%      25%      50%      75%    97.5%  Rhat n.eff
# AFT[2]       0.780   0.131    0.549    0.687    0.771    0.865    1.061 1.004   780
# AFT[3]       0.983   0.124    0.763    0.897    0.975    1.058    1.249 1.001  4000
# AFT[4]       0.898   0.110    0.701    0.822    0.892    0.966    1.131 1.002  3500
# AFT[5]       0.918   0.091    0.759    0.855    0.912    0.970    1.117 1.004  1900
# HR[2]        0.732   0.157    0.469    0.618    0.718    0.832    1.076 1.002  2000
# HR[3]        0.979   0.154    0.709    0.868    0.967    1.075    1.312 1.001  4000
# HR[4]        0.873   0.135    0.635    0.778    0.865    0.958    1.163 1.001  3200
# HR[5]        0.897   0.112    0.707    0.818    0.889    0.962    1.146 1.003  2200
# beta[1]      6.436   0.264    6.054    6.282    6.405    6.545    7.012 1.058    69
# beta[2]      0.263   0.168   -0.059    0.145    0.261    0.376    0.601 1.004   780
# beta[3]      0.025   0.125   -0.222   -0.056    0.026    0.109    0.270 1.001  4000
# beta[4]      0.115   0.122   -0.123    0.035    0.114    0.196    0.356 1.002  3500
# beta[5]      0.090   0.098   -0.110    0.030    0.092    0.156    0.276 1.004  1900
# p.crit[2]    0.947   0.224    0.000    1.000    1.000    1.000    1.000 1.010   960
# p.crit[3]    0.581   0.494    0.000    0.000    1.000    1.000    1.000 1.001  4000
# p.crit[4]    0.831   0.375    0.000    1.000    1.000    1.000    1.000 1.003  1100
# p.crit[5]    0.830   0.375    0.000    1.000    1.000    1.000    1.000 1.002  1500
# shape        1.274   0.074    1.119    1.235    1.278    1.320    1.400 1.062    56
# deviance  7755.845  47.839 7664.020 7723.739 7754.564 7785.871 7854.578 1.006   510
# 
# For each parameter, n.eff is a crude measure of effective sample size,
# and Rhat is the potential scale reduction factor (at convergence, Rhat=1).
# 
# DIC info (using the rule, pD = var(deviance)/2)
# pD = 1138.3 and DIC = 8894.2
# DIC is an estimate of expected predictive error (lower deviance is better).
# > summary(mcmcplots::as.mcmc.rjags(M0FEMHemoStroke_20200704))
# 
# Iterations = 1:97389
# Thinning interval = 97 
# Number of chains = 4 
# Sample size per chain = 1005 
# 
# 1. Empirical mean and standard deviation for each variable,
# plus standard error of the mean:
#   
#               Mean       SD Naive SE Time-series SE
# AFT[2]    7.799e-01  0.13116 0.002069       0.004767
# AFT[3]    9.828e-01  0.12436 0.001961       0.005866
# AFT[4]    8.978e-01  0.10953 0.001727       0.004883
# AFT[5]    9.183e-01  0.09082 0.001432       0.005201
# HR[2]     7.323e-01  0.15664 0.002470       0.005832
# HR[3]     9.785e-01  0.15432 0.002434       0.007165
# HR[4]     8.728e-01  0.13506 0.002130       0.006134
# HR[5]     8.972e-01  0.11181 0.001763       0.006485
# beta[1]   6.436e+00  0.26383 0.004161       0.037299
# beta[2]   2.627e-01  0.16837 0.002656       0.006170
# beta[3]   2.521e-02  0.12520 0.001975       0.005953
# beta[4]   1.152e-01  0.12157 0.001917       0.005301
# beta[5]   9.005e-02  0.09757 0.001539       0.005426
# deviance  7.756e+03 47.83928 0.754521       1.889654
# p.crit[2] 9.470e-01  0.22403 0.003533       0.005183
# p.crit[3] 5.806e-01  0.49352 0.007784       0.018633
# p.crit[4] 8.311e-01  0.37472 0.005910       0.013632
# p.crit[5] 8.303e-01  0.37537 0.005920       0.018546
# shape     1.274e+00  0.07370 0.001162       0.009604
# 
# 2. Quantiles for each variable:
#   
#   2.5%        25%       50%       75%     97.5%
# AFT[2]       0.54853    0.68673 7.705e-01    0.8651    1.0613
# AFT[3]       0.76333    0.89686 9.747e-01    1.0581    1.2487
# AFT[4]       0.70070    0.82219 8.923e-01    0.9660    1.1308
# AFT[5]       0.75882    0.85520 9.118e-01    0.9702    1.1168
# HR[2]        0.46859    0.61828 7.177e-01    0.8322    1.0756
# HR[3]        0.70872    0.86836 9.675e-01    1.0750    1.3120
# HR[4]        0.63498    0.77811 8.647e-01    0.9575    1.1635
# HR[5]        0.70653    0.81782 8.887e-01    0.9622    1.1463
# beta[1]      6.05363    6.28213 6.405e+00    6.5450    7.0117
# beta[2]     -0.05945    0.14486 2.607e-01    0.3758    0.6005
# beta[3]     -0.22207   -0.05646 2.568e-02    0.1089    0.2701
# beta[4]     -0.12294    0.03461 1.139e-01    0.1958    0.3557
# beta[5]     -0.11048    0.03022 9.228e-02    0.1564    0.2760
# deviance  7664.02013 7723.73878 7.755e+03 7785.8707 7854.5775
# p.crit[2]    0.00000    1.00000 1.000e+00    1.0000    1.0000
# p.crit[3]    0.00000    0.00000 1.000e+00    1.0000    1.0000
# p.crit[4]    0.00000    1.00000 1.000e+00    1.0000    1.0000
# p.crit[5]    0.00000    1.00000 1.000e+00    1.0000    1.0000
# shape        1.11855    1.23498 1.278e+00    1.3203    1.4000

#######################Model 1 age-adj######################
#                                                          #
#                      Model 1 age-adj                     #
#                                                          #
#**********************************************************#

Age <- MData_fem$Age


Agehemo.weibull.model <- function() {
  for(i in 1:54999){
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

MILKdataFEM <- c("t",  "c", "Mlkfre", "is.censored", "Age")
age.params <- c("AFT", "HR", "p.crit")
start.time <- Sys.time()
jagsfit <- jags.parallel(data=MILKdataFEM,  parameters.to.save = age.params, 
                         n.iter=100000, n.burnin=(5000/2), n.chains = 4,
                         model.file=Agehemo.weibull.model)
end.time <- Sys.time()
end.time - start.time #Time difference of 1.255446 days

M1FEMHemoStroke_20200704 <- jagsfit
print(M1FEMHemoStroke_20200704)
summary(mcmcplots::as.mcmc.rjags(M1FEMHemoStroke_20200704))
# save.image(file = "data/JACCmilkstrokewithHemo.Rdata")
save(M0FEMHemoStroke_20200704, 
     M1FEMHemoStroke_20200704, 
     file = "HemoStrokeFEM.RData")

# Inference for Bugs model at "Agehemo.weibull.model", fit using jags,
# 4 chains, each with 1e+05 iterations (first 2500 discarded), n.thin = 97
# n.sims = 4020 iterations saved
#            mu.vect sd.vect     2.5%      25%      50%      75%    97.5%  Rhat n.eff
# AFT[2]       0.879   0.138    0.634    0.781    0.870    0.964    1.172 1.001  3600
# AFT[3]       1.120   0.132    0.897    1.032    1.108    1.199    1.411 1.005   740
# AFT[4]       1.043   0.128    0.822    0.953    1.035    1.122    1.319 1.001  2500
# AFT[5]       0.949   0.091    0.800    0.885    0.941    1.001    1.141 1.003  2100
# HR[2]        0.842   0.182    0.538    0.710    0.825    0.952    1.243 1.001  3200
# HR[3]        1.171   0.184    0.858    1.045    1.153    1.282    1.575 1.003  1100
# HR[4]        1.063   0.177    0.764    0.935    1.049    1.172    1.446 1.001  4000
# HR[5]        0.931   0.121    0.732    0.844    0.920    1.002    1.199 1.002  1900
# p.crit[2]    0.816   0.387    0.000    1.000    1.000    1.000    1.000 1.001  3400
# p.crit[3]    0.169   0.375    0.000    0.000    0.000    0.000    1.000 1.002  1400
# p.crit[4]    0.389   0.487    0.000    0.000    0.000    1.000    1.000 1.001  4000
# p.crit[5]    0.746   0.435    0.000    0.000    1.000    1.000    1.000 1.002  2500
# deviance  7423.580  46.934 7334.193 7391.375 7423.324 7454.439 7517.134 1.007   440
# 
# For each parameter, n.eff is a crude measure of effective sample size,
# and Rhat is the potential scale reduction factor (at convergence, Rhat=1).
# 
# DIC info (using the rule, pD = var(deviance)/2)
# pD = 1094.6 and DIC = 8518.2
# DIC is an estimate of expected predictive error (lower deviance is better).
# 
# Iterations = 1:97389
# Thinning interval = 97 
# Number of chains = 4 
# Sample size per chain = 1005 
# 
# 1. Empirical mean and standard deviation for each variable,
# plus standard error of the mean:
#   
#                Mean       SD Naive SE Time-series SE
# AFT[2]       0.8792  0.13795 0.002176       0.004909
# AFT[3]       1.1198  0.13216 0.002084       0.006786
# AFT[4]       1.0434  0.12836 0.002025       0.006285
# AFT[5]       0.9486  0.09063 0.001429       0.005239
# HR[2]        0.8419  0.18195 0.002870       0.006560
# HR[3]        1.1707  0.18411 0.002904       0.008940
# HR[4]        1.0626  0.17680 0.002789       0.008474
# HR[5]        0.9306  0.12108 0.001910       0.006958
# deviance  7423.5796 46.93423 0.740247       1.328401
# p.crit[2]    0.8162  0.38739 0.006110       0.011702
# p.crit[3]    0.1694  0.37515 0.005917       0.013759
# p.crit[4]    0.3886  0.48748 0.007689       0.019956
# p.crit[5]    0.7463  0.43520 0.006864       0.021932
# 
# 2. Quantiles for each variable:
#   
#                2.5%       25%       50%       75%    97.5%
# AFT[2]       0.6339    0.7811    0.8695    0.9644    1.172
# AFT[3]       0.8966    1.0316    1.1081    1.1986    1.411
# AFT[4]       0.8223    0.9526    1.0348    1.1218    1.319
# AFT[5]       0.7997    0.8847    0.9412    1.0011    1.141
# HR[2]        0.5379    0.7096    0.8246    0.9517    1.243
# HR[3]        0.8583    1.0445    1.1528    1.2821    1.575
# HR[4]        0.7638    0.9350    1.0487    1.1718    1.446
# HR[5]        0.7319    0.8442    0.9197    1.0016    1.199
# deviance  7334.1930 7391.3748 7423.3238 7454.4386 7517.134
# p.crit[2]    0.0000    1.0000    1.0000    1.0000    1.000
# p.crit[3]    0.0000    0.0000    0.0000    0.0000    1.000
# p.crit[4]    0.0000    0.0000    0.0000    1.0000    1.000
# p.crit[5]    0.0000    0.0000    1.0000    1.0000    1.000

#######################Model 2 Mult-ad######################
#                                                          #
#                      Model 2 Mult-ad                     #
#                                                          #
#**********************************************************#



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




M2hemo.weibull.model <- function() {
  for(i in 1:54999){
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

MILKdataFEM <- c("t",  "c", "Mlkfre", "is.censored", "Age", "Smoking", 
                 "Alc_Fre", "BMIgrp"#, 
                 # "DM_hist", "HT_hist", "KID_hist", 
                 # "LIV_hist", "Exercise", "Slepgrp", "Spi", "Cofe", 
                 # "Educgrp", "Gretea", "Fru"
)
M2.params <- c("AFT", "HR", "p.crit")
start.time <- Sys.time()
jagsfit <- jags.parallel(data=MILKdataFEM,  parameters.to.save = M2.params, 
                         n.iter=100000, n.burnin=(3000/2), n.chains = 4,
                         model.file=M2hemo.weibull.model)
end.time <- Sys.time()
end.time - start.time #Time difference of 3.507339 days

M2FEMHemoStroke_20200704 <- jagsfit
print(M2FEMHemoStroke_20200704)
summary(mcmcplots::as.mcmc.rjags(M2FEMHemoStroke_20200704))
# save.image(file = "data/JACCmilkstrokewithHemo.Rdata")

save(M0FEMHemoStroke_20200704, 
     M1FEMHemoStroke_20200704,
     M2FEMHemoStroke_20200704, file = "IscheStrokeMen.RData")

# Inference for Bugs model at "M2hemo.weibull.model", fit using jags,
# 4 chains, each with 1e+05 iterations (first 1500 discarded), n.thin = 98
# n.sims = 4020 iterations saved
#            mu.vect sd.vect     2.5%      25%      50%      75%    97.5%  Rhat n.eff
# AFT[2]       0.931   0.236    0.642    0.807    0.902    1.010    1.335 1.030   360
# AFT[3]       1.226   0.385    0.926    1.073    1.160    1.264    1.983 1.061   710
# AFT[4]       1.138   0.330    0.868    1.001    1.081    1.180    1.826 1.072   530
# AFT[5]       1.035   0.248    0.831    0.935    0.993    1.059    1.547 1.080   390
# HR[2]        0.899   0.223    0.553    0.747    0.871    1.013    1.396 1.008   560
# HR[3]        1.264   0.260    0.902    1.101    1.223    1.370    1.903 1.009  1100
# HR[4]        1.150   0.234    0.827    1.001    1.112    1.251    1.741 1.016   560
# HR[5]        1.021   0.186    0.779    0.913    0.991    1.080    1.509 1.021   560
# p.crit[2]    0.732   0.443    0.000    0.000    1.000    1.000    1.000 1.003  1200
# p.crit[3]    0.095   0.294    0.000    0.000    0.000    0.000    1.000 1.005  1100
# p.crit[4]    0.248   0.432    0.000    0.000    0.000    0.000    1.000 1.005   590
# p.crit[5]    0.533   0.499    0.000    0.000    1.000    1.000    1.000 1.008   350
# deviance  7402.158  78.744 7300.863 7358.306 7391.527 7426.013 7605.092 1.045   580
# 
# For each parameter, n.eff is a crude measure of effective sample size,
# and Rhat is the potential scale reduction factor (at convergence, Rhat=1).
# 
# DIC info (using the rule, pD = var(deviance)/2)
# pD = 3086.2 and DIC = 10488.4
# DIC is an estimate of expected predictive error (lower deviance is better).
# > summary(mcmcplots::as.mcmc.rjags(M2FEMHemoStroke_20200704))
# 
# Iterations = 1:98393
# Thinning interval = 98 
# Number of chains = 4 
# Sample size per chain = 1005 
# 
# 1. Empirical mean and standard deviation for each variable,
# plus standard error of the mean:
#   
#                Mean      SD Naive SE Time-series SE
# AFT[2]    9.306e-01  0.2357 0.003718       0.018172
# AFT[3]    1.226e+00  0.3851 0.006073       0.039604
# AFT[4]    1.138e+00  0.3301 0.005206       0.030179
# AFT[5]    1.035e+00  0.2480 0.003912       0.025801
# HR[2]     8.986e-01  0.2235 0.003525       0.015941
# HR[3]     1.264e+00  0.2601 0.004102       0.026417
# HR[4]     1.150e+00  0.2343 0.003695       0.026560
# HR[5]     1.021e+00  0.1855 0.002926       0.021356
# deviance  7.402e+03 78.7444 1.241958       8.643872
# p.crit[2] 7.321e-01  0.4429 0.006986       0.022380
# p.crit[3] 9.527e-02  0.2936 0.004631       0.008107
# p.crit[4] 2.478e-01  0.4318 0.006810       0.012885
# p.crit[5] 5.331e-01  0.4990 0.007870       0.022618
# 
# 2. Quantiles for each variable:
#   
#                2.5%       25%       50%      75%    97.5%
# AFT[2]       0.6417    0.8072    0.9021    1.010    1.335
# AFT[3]       0.9262    1.0731    1.1598    1.264    1.983
# AFT[4]       0.8682    1.0010    1.0812    1.180    1.826
# AFT[5]       0.8311    0.9346    0.9932    1.059    1.547
# HR[2]        0.5534    0.7471    0.8708    1.013    1.396
# HR[3]        0.9022    1.1008    1.2230    1.370    1.903
# HR[4]        0.8270    1.0014    1.1118    1.251    1.741
# HR[5]        0.7787    0.9133    0.9907    1.080    1.509
# deviance  7300.8629 7358.3064 7391.5266 7426.013 7605.092
# p.crit[2]    0.0000    0.0000    1.0000    1.000    1.000
# p.crit[3]    0.0000    0.0000    0.0000    0.000    1.000
# p.crit[4]    0.0000    0.0000    0.0000    0.000    1.000
# p.crit[5]    0.0000    0.0000    1.0000    1.000    1.000

