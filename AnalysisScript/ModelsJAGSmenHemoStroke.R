# HemoStroke == "I60_2"

library(tidyverse)

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
is.censored <- MData_men$HemoStroke != "I60_2"

# define followup time for events
t <- if_else(!is.censored, MData_men$followpy, 0)
t <- na_if(t, 0)
t


# define followup time for censored 
c <- if_else(is.censored, MData_men$followpy, 0)

c[!is.censored] <- t[!is.censored] + 1




jacchemo.weibull.model0 <- function() {
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
                         model.file=jacchemo.weibull.model0)
end.time <- Sys.time()
end.time - start.time #Time difference of 9.715389 hours
M0menHemoStroke_20200629 <- jagsfit
print(M0menHemoStroke_20200629)
summary(mcmcplots::as.mcmc.rjags(M0menHemoStroke_20200629))
save.image(file = "data/JACCmilkstrokewithHemo.Rdata")
# save.image(file = "data/JACCmilkstroke.Rdata")

# Inference for Bugs model at "jacchemo.weibull.model0", fit using jags,
# 4 chains, each with 1e+05 iterations (first 2500 discarded), n.thin = 97
# n.sims = 4020 iterations saved
# mu.vect sd.vect     2.5%      25%      50%      75%    97.5%  Rhat n.eff
# AFT[2]       1.028   0.165    0.735    0.914    1.012    1.129    1.378 1.004   790
# AFT[3]       0.851   0.124    0.632    0.764    0.843    0.929    1.116 1.010   260
# AFT[4]       0.865   0.126    0.645    0.778    0.856    0.942    1.140 1.004   710
# AFT[5]       0.975   0.110    0.784    0.899    0.963    1.042    1.223 1.008   330
# HR[2]        1.034   0.196    0.697    0.899    1.015    1.154    1.460 1.004   770
# HR[3]        0.827   0.143    0.579    0.726    0.815    0.915    1.135 1.009   290
# HR[4]        0.843   0.146    0.597    0.741    0.830    0.931    1.168 1.004   790
# HR[5]        0.970   0.129    0.748    0.882    0.956    1.050    1.259 1.008   340
# beta[1]      6.545   0.229    6.165    6.392    6.531    6.661    7.030 1.021   360
# beta[2]     -0.015   0.159   -0.320   -0.121   -0.012    0.090    0.308 1.004   790
# beta[3]      0.172   0.145   -0.110    0.074    0.171    0.269    0.459 1.010   260
# beta[4]      0.155   0.144   -0.131    0.060    0.155    0.250    0.439 1.004   710
# beta[5]      0.032   0.111   -0.201   -0.041    0.037    0.106    0.244 1.008   330
# p.crit[2]    0.472   0.499    0.000    0.000    0.000    1.000    1.000 1.002  1300
# p.crit[3]    0.884   0.320    0.000    1.000    1.000    1.000    1.000 1.006   900
# p.crit[4]    0.863   0.344    0.000    1.000    1.000    1.000    1.000 1.004  1200
# p.crit[5]    0.631   0.483    0.000    0.000    1.000    1.000    1.000 1.004   730
# shape        1.190   0.061    1.065    1.152    1.190    1.229    1.304 1.025   210
# deviance  6293.186  42.019 6212.679 6264.635 6291.913 6321.312 6378.146 1.003  1200
# 
# For each parameter, n.eff is a crude measure of effective sample size,
# and Rhat is the potential scale reduction factor (at convergence, Rhat=1).
# 
# DIC info (using the rule, pD = var(deviance)/2)
# pD = 881.2 and DIC = 7174.4
# DIC is an estimate of expected predictive error (lower deviance is better).

# Iterations = 1:97389
# Thinning interval = 97 
# Number of chains = 4 
# Sample size per chain = 1005 
# 
# 1. Empirical mean and standard deviation for each variable,
# plus standard error of the mean:
#   
#                 Mean       SD Naive SE Time-series SE
# AFT[2]       1.02769  0.16489 0.002601       0.005675
# AFT[3]       0.85116  0.12375 0.001952       0.004883
# AFT[4]       0.86513  0.12594 0.001986       0.004930
# AFT[5]       0.97450  0.10995 0.001734       0.005410
# HR[2]        1.03447  0.19587 0.003089       0.006645
# HR[3]        0.82729  0.14254 0.002248       0.005667
# HR[4]        0.84327  0.14569 0.002298       0.005762
# HR[5]        0.97003  0.12921 0.002038       0.006230
# beta[1]      6.54491  0.22874 0.003608       0.025734
# beta[2]     -0.01463  0.15930 0.002513       0.005468
# beta[3]      0.17166  0.14506 0.002288       0.005740
# beta[4]      0.15531  0.14433 0.002276       0.005584
# beta[5]      0.03207  0.11140 0.001757       0.005452
# deviance  6293.18638 42.01875 0.662720       1.332048
# p.crit[2]    0.47189  0.49927 0.007875       0.015085
# p.crit[3]    0.88408  0.32017 0.005050       0.009483
# p.crit[4]    0.86318  0.34370 0.005421       0.010291
# p.crit[5]    0.63109  0.48257 0.007611       0.020265
# shape        1.18970  0.06087 0.000960       0.006748
# 
# 2. Quantiles for each variable:
#   
#                2.5%        25%        50%       75%     97.5%
# AFT[2]       0.7349    0.91415    1.01242 1.129e+00    1.3776
# AFT[3]       0.6321    0.76412    0.84318 9.287e-01    1.1160
# AFT[4]       0.6450    0.77848    0.85599 9.417e-01    1.1401
# AFT[5]       0.7838    0.89935    0.96342 1.042e+00    1.2231
# HR[2]        0.6966    0.89860    1.01482 1.154e+00    1.4601
# HR[3]        0.5795    0.72611    0.81539 9.151e-01    1.1351
# HR[4]        0.5968    0.74098    0.83040 9.309e-01    1.1682
# HR[5]        0.7483    0.88195    0.95637 1.050e+00    1.2588
# beta[1]      6.1649    6.39154    6.53084 6.661e+00    7.0303
# beta[2]     -0.3203   -0.12104   -0.01234 8.976e-02    0.3080
# beta[3]     -0.1097    0.07400    0.17058 2.690e-01    0.4587
# beta[4]     -0.1312    0.06010    0.15550 2.504e-01    0.4386
# beta[5]     -0.2014   -0.04103    0.03726 1.061e-01    0.2436
# deviance  6212.6792 6264.63497 6291.91276 6.321e+03 6378.1464
# p.crit[2]    0.0000    0.00000    0.00000 1.000e+00    1.0000
# p.crit[3]    0.0000    1.00000    1.00000 1.000e+00    1.0000
# p.crit[4]    0.0000    1.00000    1.00000 1.000e+00    1.0000
# p.crit[5]    0.0000    0.00000    1.00000 1.000e+00    1.0000
# shape        1.0653    1.15227    1.18998 1.229e+00    1.3038





#######################Model 1 age-adj######################
#                                                          #
#                      Model 1 age-adj                     #
#                                                          #
#**********************************************************#

Age <- MData_men$Age


Agehemo.weibull.model <- function() {
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
                         model.file=Agehemo.weibull.model)
end.time <- Sys.time()
end.time - start.time #Time difference of 13.87173 hours

M1menHemoStroke_20200630 <- jagsfit
print(M1menHemoStroke_20200630)
summary(mcmcplots::as.mcmc.rjags(M1menHemoStroke_20200630))
save.image(file = "data/JACCmilkstrokewithHemo.Rdata")

# Inference for Bugs model at "Agehemo.weibull.model", fit using jags,
# 4 chains, each with 1e+05 iterations (first 2500 discarded), n.thin = 97
# n.sims = 4020 iterations saved
#            mu.vect sd.vect     2.5%      25%      50%      75%    97.5%  Rhat n.eff
# AFT[2]       1.087   0.171    0.796    0.972    1.072    1.183    1.455 1.002  2800
# AFT[3]       0.910   0.133    0.699    0.823    0.898    0.978    1.200 1.002  4000
# AFT[4]       0.923   0.132    0.705    0.835    0.915    0.994    1.195 1.003  2200
# AFT[5]       0.904   0.103    0.741    0.841    0.896    0.955    1.114 1.009  1300
# HR[2]        1.113   0.210    0.745    0.964    1.095    1.241    1.586 1.001  4000
# HR[3]        0.884   0.158    0.628    0.774    0.869    0.972    1.252 1.001  4000
# HR[4]        0.900   0.157    0.634    0.789    0.890    0.992    1.244 1.001  4000
# HR[5]        0.875   0.119    0.674    0.796    0.866    0.942    1.145 1.002  2300
# p.crit[2]    0.316   0.465    0.000    0.000    0.000    1.000    1.000 1.001  4000
# p.crit[3]    0.797   0.402    0.000    1.000    1.000    1.000    1.000 1.001  3900
# p.crit[4]    0.766   0.423    0.000    1.000    1.000    1.000    1.000 1.001  3800
# p.crit[5]    0.876   0.329    0.000    1.000    1.000    1.000    1.000 1.001  4000
# deviance  6078.934  46.010 5995.277 6049.644 6078.106 6106.763 6167.160 1.013  3600

# Iterations = 1:97389
# Thinning interval = 97 
# Number of chains = 4 
# Sample size per chain = 1005 
# 
# 1. Empirical mean and standard deviation for each variable,
# plus standard error of the mean:
#   
#                Mean      SD Naive SE Time-series SE
# AFT[2]       1.0867  0.1710 0.002698       0.009351
# AFT[3]       0.9097  0.1335 0.002105       0.006623
# AFT[4]       0.9228  0.1324 0.002088       0.007032
# AFT[5]       0.9036  0.1026 0.001618       0.007499
# HR[2]        1.1125  0.2104 0.003318       0.009051
# HR[3]        0.8840  0.1581 0.002494       0.006712
# HR[4]        0.9005  0.1571 0.002477       0.007352
# HR[5]        0.8753  0.1190 0.001877       0.007663
# deviance  6078.9344 46.0101 0.725672       2.181546
# p.crit[2]    0.3162  0.4650 0.007335       0.012503
# p.crit[3]    0.7968  0.4025 0.006348       0.018283
# p.crit[4]    0.7659  0.4235 0.006679       0.015397
# p.crit[5]    0.8761  0.3295 0.005197       0.017055
# 
# 2. Quantiles for each variable:
#   
#                2.5%       25%       50%       75%    97.5%
# AFT[2]       0.7958    0.9724    1.0718    1.1826    1.455
# AFT[3]       0.6994    0.8226    0.8980    0.9785    1.200
# AFT[4]       0.7048    0.8348    0.9145    0.9939    1.195
# AFT[5]       0.7414    0.8409    0.8962    0.9552    1.114
# HR[2]        0.7453    0.9639    1.0951    1.2414    1.586
# HR[3]        0.6281    0.7736    0.8690    0.9718    1.252
# HR[4]        0.6337    0.7893    0.8902    0.9918    1.244
# HR[5]        0.6736    0.7957    0.8659    0.9424    1.145
# deviance  5995.2772 6049.6438 6078.1063 6106.7635 6167.160
# p.crit[2]    0.0000    0.0000    0.0000    1.0000    1.000
# p.crit[3]    0.0000    1.0000    1.0000    1.0000    1.000
# p.crit[4]    0.0000    1.0000    1.0000    1.0000    1.000
# p.crit[5]    0.0000    1.0000    1.0000    1.0000    1.000



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
end.time - start.time  #Time difference of 2.548741 days
M1menHemoStroke_20200701 <- jagsfit
print(M1menHemoStroke_20200701)
summary(mcmcplots::as.mcmc.rjags(M1menHemoStroke_20200701))
save.image(file = "data/JACCmilkstrokewithHemo.Rdata")


# Iterations = 1:98393
# Thinning interval = 98 
# Number of chains = 4 
# Sample size per chain = 1005 
# 
# 1. Empirical mean and standard deviation for each variable,
# plus standard error of the mean:
#   
#   Mean      SD Naive SE Time-series SE
# AFT[2]       1.1111  0.1985 0.003130       0.010852
# AFT[3]       0.9333  0.1467 0.002314       0.007357
# AFT[4]       0.9628  0.1592 0.002511       0.008552
# AFT[5]       0.9596  0.1311 0.002067       0.007722
# HR[2]        1.1373  0.2244 0.003539       0.009477
# HR[3]        0.9169  0.1654 0.002609       0.007575
# HR[4]        0.9527  0.1823 0.002875       0.009809
# HR[5]        0.9477  0.1428 0.002252       0.008262
# deviance  6060.8100 59.9484 0.945507       5.225167
# p.crit[2]    0.2881  0.4529 0.007143       0.012517
# p.crit[3]    0.7241  0.4470 0.007050       0.016330
# p.crit[4]    0.6435  0.4790 0.007555       0.020082
# p.crit[5]    0.6930  0.4613 0.007276       0.020135
# 
# 2. Quantiles for each variable:
#   
#               2.5%       25%       50%      75%    97.5%
# AFT[2]       0.7950    0.9848    1.0928    1.214    1.575
# AFT[3]       0.6957    0.8419    0.9210    1.011    1.249
# AFT[4]       0.7085    0.8605    0.9498    1.043    1.340
# AFT[5]       0.7592    0.8841    0.9489    1.021    1.247
# HR[2]        0.7524    0.9815    1.1202    1.272    1.613
# HR[3]        0.6294    0.8049    0.9014    1.013    1.287
# HR[4]        0.6531    0.8248    0.9374    1.053    1.372
# HR[5]        0.7104    0.8540    0.9362    1.027    1.276
# deviance  5971.2341 6025.2626 6053.7141 6085.181 6206.865
# p.crit[2]    0.0000    0.0000    0.0000    1.000    1.000
# p.crit[3]    0.0000    0.0000    1.0000    1.000    1.000
# p.crit[4]    0.0000    0.0000    1.0000    1.000    1.000
# p.crit[5]    0.0000    0.0000    1.0000    1.000    1.000


# Inference for Bugs model at "M2hemo.weibull.model", fit using jags,
# 4 chains, each with 1e+05 iterations (first 1500 discarded), n.thin = 98
# n.sims = 4020 iterations saved
#            mu.vect sd.vect     2.5%      25%      50%      75%    97.5%  Rhat n.eff
# AFT[2]       1.111   0.198    0.795    0.985    1.093    1.214    1.575 1.002  1600
# AFT[3]       0.933   0.147    0.696    0.842    0.921    1.011    1.249 1.005   590
# AFT[4]       0.963   0.159    0.709    0.861    0.950    1.043    1.340 1.003   940
# AFT[5]       0.960   0.131    0.759    0.884    0.949    1.021    1.247 1.007   660
# HR[2]        1.137   0.224    0.752    0.982    1.120    1.272    1.613 1.002  1400
# HR[3]        0.917   0.165    0.629    0.805    0.901    1.013    1.287 1.005   530
# HR[4]        0.953   0.182    0.653    0.825    0.937    1.053    1.372 1.004   770
# HR[5]        0.948   0.143    0.710    0.854    0.936    1.027    1.276 1.006   500
# p.crit[2]    0.288   0.453    0.000    0.000    0.000    1.000    1.000 1.001  3600
# p.crit[3]    0.724   0.447    0.000    0.000    1.000    1.000    1.000 1.003  1000
# p.crit[4]    0.644   0.479    0.000    0.000    1.000    1.000    1.000 1.002  2000
# p.crit[5]    0.693   0.461    0.000    0.000    1.000    1.000    1.000 1.005   570
# deviance  6060.810  59.948 5971.234 6025.263 6053.714 6085.181 6206.865 1.002  4000
# 
# For each parameter, n.eff is a crude measure of effective sample size,
# and Rhat is the potential scale reduction factor (at convergence, Rhat=1).
# 
# DIC info (using the rule, pD = var(deviance)/2)
# pD = 1797.8 and DIC = 7858.6
# DIC is an estimate of expected predictive error (lower deviance is better).