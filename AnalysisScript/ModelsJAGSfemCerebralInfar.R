# IschemStroke == "I63"

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
is.censored <- MData_fem$IscheStroke != "I63"

# define followup time for events
t <- if_else(!is.censored, MData_fem$followpy, 0)
t <- na_if(t, 0)
t


# define followup time for censored 
c <- if_else(is.censored, MData_fem$followpy, 0)

c[!is.censored] <- t[!is.censored] + 1




jaccIschemic.weibull.model0 <- function() {
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

MILKdatafem <- c("t",  "c", "Mlkfre", "is.censored")
m0.params <- c("beta", "AFT", "HR", "p.crit", "shape")
library(R2jags)


start.time <- Sys.time()
jagsfit <- jags.parallel(data=MILKdatafem,  parameters.to.save = m0.params, 
                         n.iter=100000, n.burnin=(5000/2), n.chains = 4, #n.thin = 100,
                         model.file=jaccIschemic.weibull.model0)
end.time <- Sys.time()
end.time - start.time #
M0femIscheStroke_20200707 <- jagsfit
print(M0femIscheStroke_20200707)
summary(mcmcplots::as.mcmc.rjags(M0femIscheStroke_20200707))

save(M0femIscheStroke_20200707, file = "IscheStrokefem.RData")

# Inference for Bugs model at "jaccIschemic.weibull.model0", fit using jags,
# 4 chains, each with 1e+05 iterations (first 2500 discarded), n.thin = 97
# n.sims = 4020 iterations saved
#             mu.vect sd.vect     2.5%      25%      50%      75%    97.5%  Rhat n.eff
# AFT[2]       1.013   0.125    0.789    0.927    1.007    1.090    1.267 1.001  3700
# AFT[3]       0.902   0.089    0.746    0.842    0.896    0.957    1.099 1.007   480
# AFT[4]       0.745   0.078    0.600    0.690    0.742    0.794    0.906 1.007   420
# AFT[5]       0.863   0.064    0.746    0.817    0.862    0.904    0.998 1.006   540
# HR[2]        1.028   0.207    0.671    0.881    1.013    1.156    1.468 1.001  4000
# HR[3]        0.845   0.138    0.613    0.749    0.830    0.929    1.159 1.006   580
# HR[4]        0.614   0.108    0.428    0.538    0.604    0.679    0.850 1.006   500
# HR[5]        0.783   0.097    0.615    0.713    0.779    0.843    0.996 1.005   760
# beta[1]      5.639   0.211    5.361    5.523    5.609    5.712    6.085 1.075   160
# beta[2]     -0.005   0.122   -0.236   -0.086   -0.007    0.075    0.237 1.001  3700
# beta[3]      0.108   0.097   -0.094    0.044    0.110    0.172    0.293 1.007   480
# beta[4]      0.300   0.105    0.098    0.230    0.298    0.371    0.511 1.007   420
# beta[5]      0.150   0.074    0.002    0.101    0.149    0.202    0.292 1.006   540
# p.crit[2]    0.471   0.499    0.000    0.000    0.000    1.000    1.000 1.001  3500
# p.crit[3]    0.864   0.343    0.000    1.000    1.000    1.000    1.000 1.008   570
# p.crit[4]    0.997   0.052    1.000    1.000    1.000    1.000    1.000 1.137  1100
# p.crit[5]    0.976   0.153    1.000    1.000    1.000    1.000    1.000 1.014  1500
# shape        1.675   0.100    1.463    1.624    1.683    1.736    1.838 1.060   140
# deviance  6570.107  44.278 6487.846 6539.675 6568.745 6598.676 6658.043 1.003  1800
# 
# For each parameter, n.eff is a crude measure of effective sample size,
# and Rhat is the potential scale reduction factor (at convergence, Rhat=1).
# 
# DIC info (using the rule, pD = var(deviance)/2)
# pD = 979.4 and DIC = 7549.5
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
#                 Mean       SD Naive SE Time-series SE
# AFT[2]     1.013e+00  0.12480 0.001968       0.004497
# AFT[3]     9.022e-01  0.08863 0.001398       0.003756
# AFT[4]     7.447e-01  0.07826 0.001234       0.003468
# AFT[5]     8.627e-01  0.06397 0.001009       0.003456
# HR[2]      1.028e+00  0.20650 0.003257       0.007163
# HR[3]      8.452e-01  0.13830 0.002181       0.005846
# HR[4]      6.143e-01  0.10847 0.001711       0.005237
# HR[5]      7.834e-01  0.09679 0.001527       0.005318
# beta[1]    5.639e+00  0.21075 0.003324       0.030263
# beta[2]   -5.312e-03  0.12227 0.001928       0.004291
# beta[3]    1.077e-01  0.09726 0.001534       0.004106
# beta[4]    3.002e-01  0.10502 0.001656       0.004662
# beta[5]    1.504e-01  0.07389 0.001165       0.003997
# deviance   6.570e+03 44.27827 0.698357       1.682586
# p.crit[2]  4.711e-01  0.49923 0.007874       0.014440
# p.crit[3]  8.637e-01  0.34317 0.005412       0.012422
# p.crit[4]  9.973e-01  0.05224 0.000824       0.001539
# p.crit[5]  9.761e-01  0.15270 0.002408       0.005745
# shape      1.675e+00  0.10036 0.001583       0.014232
# 
# 2. Quantiles for each variable:
#   
#                 2.5%        25%        50%       75%     97.5%
# AFT[2]     7.889e-01    0.92744  1.007e+00 1.090e+00    1.2668
# AFT[3]     7.461e-01    0.84184  8.956e-01 9.566e-01    1.0987
# AFT[4]     6.002e-01    0.69035  7.422e-01 7.943e-01    0.9064
# AFT[5]     7.465e-01    0.81674  8.618e-01 9.039e-01    0.9975
# HR[2]      6.709e-01    0.88094  1.013e+00 1.156e+00    1.4681
# HR[3]      6.130e-01    0.74903  8.295e-01 9.285e-01    1.1588
# HR[4]      4.282e-01    0.53822  6.044e-01 6.788e-01    0.8495
# HR[5]      6.148e-01    0.71303  7.787e-01 8.434e-01    0.9965
# beta[1]    5.361e+00    5.52297  5.609e+00 5.712e+00    6.0855
# beta[2]   -2.365e-01   -0.08646 -7.439e-03 7.533e-02    0.2372
# beta[3]   -9.416e-02    0.04437  1.103e-01 1.722e-01    0.2928
# beta[4]    9.824e-02    0.23026  2.981e-01 3.705e-01    0.5105
# beta[5]    2.477e-03    0.10099  1.487e-01 2.024e-01    0.2924
# deviance   6.488e+03 6539.67522  6.569e+03 6.599e+03 6658.0425
# p.crit[2]  0.000e+00    0.00000  0.000e+00 1.000e+00    1.0000
# p.crit[3]  0.000e+00    1.00000  1.000e+00 1.000e+00    1.0000
# p.crit[4]  1.000e+00    1.00000  1.000e+00 1.000e+00    1.0000
# p.crit[5]  1.000e+00    1.00000  1.000e+00 1.000e+00    1.0000
# shape      1.463e+00    1.62431  1.683e+00 1.736e+00    1.8384

#######################Model 1 age-adj######################
#                                                          #
#                      Model 1 age-adj                     #
#                                                          #
#**********************************************************#

Age <- MData_fem$Age


AgeIsche.weibull.model <- function() {
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

MILKdatafem <- c("t",  "c", "Mlkfre", "is.censored", "Age")
age.params <- c("AFT", "HR", "p.crit")
start.time <- Sys.time()
jagsfit <- jags.parallel(data=MILKdatafem,  parameters.to.save = age.params, 
                         n.iter=100000, n.burnin=(5000/2), n.chains = 4,
                         model.file=AgeIsche.weibull.model)
end.time <- Sys.time()
end.time - start.time # Time difference of 1.136057 days

M1femIscheStroke_20200707 <- jagsfit
print(M1femIscheStroke_20200707)
summary(mcmcplots::as.mcmc.rjags(M1femIscheStroke_20200707))
save(M0femIscheStroke_20200707, 
     M1femIscheStroke_20200707, file = "IscheStrokefem.RData")

# Inference for Bugs model at "AgeIsche.weibull.model", fit using jags,
# 4 chains, each with 1e+05 iterations (first 2500 discarded), n.thin = 97
# n.sims = 4020 iterations saved
           # mu.vect sd.vect     2.5%      25%      50%      75%    97.5%  Rhat n.eff
# AFT[2]       1.214   0.321    0.945    1.074    1.151    1.241    2.077 1.099    97
# AFT[3]       1.155   0.302    0.930    1.034    1.092    1.171    1.945 1.119    78
# AFT[4]       0.982   0.185    0.792    0.890    0.949    1.017    1.482 1.100    84
# AFT[5]       0.968   0.143    0.835    0.898    0.939    0.987    1.425 1.091    85
# HR[2]        1.373   0.326    0.894    1.154    1.322    1.524    2.176 1.015   300
# HR[3]        1.247   0.277    0.866    1.069    1.193    1.360    1.954 1.028   130
# HR[4]        0.938   0.217    0.632    0.794    0.902    1.034    1.519 1.028   130
# HR[5]        0.915   0.172    0.694    0.807    0.882    0.975    1.401 1.036   120
# p.crit[2]    0.085   0.279    0.000    0.000    0.000    0.000    1.000 1.001  4000
# p.crit[3]    0.142   0.349    0.000    0.000    0.000    0.000    1.000 1.003  1500
# p.crit[4]    0.701   0.458    0.000    0.000    1.000    1.000    1.000 1.011   270
# p.crit[5]    0.794   0.404    0.000    1.000    1.000    1.000    1.000 1.015   240
# deviance  5623.653 103.484 5518.329 5574.061 5603.266 5639.152 5915.564 1.141    78
# 
# For each parameter, n.eff is a crude measure of effective sample size,
# and Rhat is the potential scale reduction factor (at convergence, Rhat=1).
# 
# DIC info (using the rule, pD = var(deviance)/2)
# pD = 5147.0 and DIC = 10770.6
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
#               Mean       SD Naive SE Time-series SE
# AFT[2]    1.214e+00   0.3205 0.005055       0.055729
# AFT[3]    1.155e+00   0.3017 0.004758       0.040716
# AFT[4]    9.825e-01   0.1845 0.002911       0.028878
# AFT[5]    9.679e-01   0.1435 0.002263       0.024386
# HR[2]     1.373e+00   0.3258 0.005139       0.023856
# HR[3]     1.247e+00   0.2766 0.004362       0.029075
# HR[4]     9.382e-01   0.2171 0.003425       0.024018
# HR[5]     9.147e-01   0.1719 0.002711       0.020346
# deviance  5.624e+03 103.4841 1.632152      17.892528
# p.crit[2] 8.532e-02   0.2794 0.004407       0.006638
# p.crit[3] 1.423e-01   0.3494 0.005511       0.011925
# p.crit[4] 7.010e-01   0.4579 0.007222       0.030021
# p.crit[5] 7.943e-01   0.4043 0.006376       0.034330
# 
# 2. Quantiles for each variable:
#   
#   2.5%       25%       50%       75%    97.5%
# AFT[2]       0.9446    1.0739    1.1510    1.2408    2.077
# AFT[3]       0.9303    1.0337    1.0925    1.1709    1.945
# AFT[4]       0.7917    0.8898    0.9495    1.0166    1.482
# AFT[5]       0.8349    0.8985    0.9391    0.9872    1.425
# HR[2]        0.8940    1.1545    1.3219    1.5240    2.176
# HR[3]        0.8659    1.0692    1.1934    1.3601    1.954
# HR[4]        0.6321    0.7936    0.9021    1.0336    1.519
# HR[5]        0.6942    0.8068    0.8818    0.9749    1.401
# deviance  5518.3289 5574.0610 5603.2663 5639.1520 5915.564
# p.crit[2]    0.0000    0.0000    0.0000    0.0000    1.000
# p.crit[3]    0.0000    0.0000    0.0000    0.0000    1.000
# p.crit[4]    0.0000    0.0000    1.0000    1.0000    1.000
# p.crit[5]    0.0000    1.0000    1.0000    1.0000    1.000

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




M2Isch.weibull.model <- function() {
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

MILKdatafem <- c("t",  "c", "Mlkfre", "is.censored", "Age", "Smoking", 
                 "Alc_Fre", "BMIgrp"#, 
                 # "DM_hist", "HT_hist", "KID_hist", 
                 # "LIV_hist", "Exercise", "Slepgrp", "Spi", "Cofe", 
                 # "Educgrp", "Gretea", "Fru"
)
M2.params <- c("AFT", "HR", "p.crit")
start.time <- Sys.time()
jagsfit <- jags.parallel(data=MILKdatafem,  parameters.to.save = M2.params, 
                         n.iter=100000, n.burnin=(3000/2), n.chains = 4,
                         model.file=M2Isch.weibull.model)
end.time <- Sys.time()
end.time - start.time # Time difference of 2.856959 days


M2femIscheStroke_20200707 <- jagsfit
print(M2femIscheStroke_20200707)
summary(mcmcplots::as.mcmc.rjags(M2femIscheStroke_20200707))
save(M0femIscheStroke_20200707, 
     M1femIscheStroke_20200707, 
     M2femIscheStroke_20200707, file = "IscheStrokefem.RData")

# Inference for Bugs model at "M2Isch.weibull.model", fit using jags,
# 4 chains, each with 1e+05 iterations (first 1500 discarded), n.thin = 98
# n.sims = 4020 iterations saved
# mu.vect sd.vect     2.5%      25%      50%      75%    97.5%  Rhat n.eff
# AFT[2]       1.191   0.191    0.944    1.085    1.163    1.253    1.616 1.069   100
# AFT[3]       1.115   0.149    0.921    1.030    1.095    1.165    1.497 1.070    92
# AFT[4]       0.955   0.118    0.780    0.887    0.943    1.007    1.208 1.025   300
# AFT[5]       0.970   0.092    0.831    0.916    0.963    1.011    1.182 1.030   180
# HR[2]        1.376   0.288    0.893    1.175    1.347    1.545    2.021 1.009   310
# HR[3]        1.214   0.217    0.851    1.060    1.197    1.342    1.703 1.012   220
# HR[4]        0.912   0.182    0.615    0.789    0.890    1.012    1.335 1.010   310
# HR[5]        0.938   0.141    0.699    0.840    0.928    1.021    1.252 1.012   230
# p.crit[2]    0.073   0.260    0.000    0.000    0.000    0.000    1.000 1.004  2000
# p.crit[3]    0.156   0.363    0.000    0.000    0.000    0.000    1.000 1.003  1200
# p.crit[4]    0.728   0.445    0.000    0.000    1.000    1.000    1.000 1.007   390
# p.crit[5]    0.700   0.458    0.000    0.000    1.000    1.000    1.000 1.007   420
# deviance  5582.978  84.023 5485.273 5538.494 5568.945 5603.275 5835.761 1.121    85
# 
# For each parameter, n.eff is a crude measure of effective sample size,
# and Rhat is the potential scale reduction factor (at convergence, Rhat=1).
# 
# DIC info (using the rule, pD = var(deviance)/2)
# pD = 3404.5 and DIC = 8987.4
# DIC is an estimate of expected predictive error (lower deviance is better).
# > summary(mcmcplots::as.mcmc.rjags(M2femIscheStroke_20200707))
# 
# Iterations = 1:98393
# Thinning interval = 98 
# Number of chains = 4 
# Sample size per chain = 1005 
# 
# 1. Empirical mean and standard deviation for each variable,
# plus standard error of the mean:
#   
#   Mean       SD Naive SE Time-series SE
# AFT[2]    1.191e+00  0.19055 0.003005       0.018679
# AFT[3]    1.115e+00  0.14923 0.002354       0.018897
# AFT[4]    9.549e-01  0.11846 0.001868       0.010140
# AFT[5]    9.701e-01  0.09219 0.001454       0.009696
# HR[2]     1.376e+00  0.28818 0.004545       0.012364
# HR[3]     1.214e+00  0.21719 0.003425       0.013040
# HR[4]     9.123e-01  0.18197 0.002870       0.011905
# HR[5]     9.384e-01  0.14118 0.002227       0.009423
# deviance  5.583e+03 84.02255 1.325204      12.053794
# p.crit[2] 7.313e-02  0.26039 0.004107       0.007334
# p.crit[3] 1.565e-01  0.36334 0.005731       0.010896
# p.crit[4] 7.279e-01  0.44512 0.007020       0.021866
# p.crit[5] 6.998e-01  0.45842 0.007230       0.025345
# 
# 2. Quantiles for each variable:
#   
#   2.5%       25%       50%      75%    97.5%
# AFT[2]       0.9437    1.0847    1.1633    1.253    1.616
# AFT[3]       0.9212    1.0301    1.0946    1.165    1.497
# AFT[4]       0.7797    0.8875    0.9429    1.007    1.208
# AFT[5]       0.8309    0.9155    0.9632    1.011    1.182
# HR[2]        0.8928    1.1752    1.3470    1.545    2.021
# HR[3]        0.8509    1.0600    1.1971    1.342    1.703
# HR[4]        0.6154    0.7893    0.8896    1.012    1.335
# HR[5]        0.6991    0.8401    0.9282    1.021    1.252
# deviance  5485.2732 5538.4942 5568.9447 5603.275 5835.761
# p.crit[2]    0.0000    0.0000    0.0000    0.000    1.000
# p.crit[3]    0.0000    0.0000    0.0000    0.000    1.000
# p.crit[4]    0.0000    0.0000    1.0000    1.000    1.000
# p.crit[5]    0.0000    0.0000    1.0000    1.000    1.000