# IschemStroke == "I63"

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
is.censored <- MData_men$IscheStroke != "I63"

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

# Inference for Bugs model at "jaccIschemic.weibull.model0", fit using jags,
# 4 chains, each with 1e+05 iterations (first 2500 discarded), n.thin = 97
# n.sims = 4020 iterations saved
#            mu.vect sd.vect     2.5%      25%      50%      75%    97.5%  Rhat n.eff
# AFT[2]       0.755   0.091    0.593    0.693    0.749    0.813    0.944 1.004   750
# AFT[3]       0.706   0.072    0.575    0.654    0.703    0.752    0.855 1.002  4000
# AFT[4]       0.744   0.075    0.609    0.692    0.741    0.790    0.896 1.002  2000
# AFT[5]       0.798   0.064    0.684    0.755    0.795    0.835    0.925 1.009   440
# HR[2]        0.657   0.121    0.450    0.575    0.646    0.731    0.919 1.004   710
# HR[3]        0.593   0.093    0.434    0.527    0.587    0.651    0.797 1.002  4000
# HR[4]        0.641   0.098    0.472    0.573    0.636    0.702    0.848 1.003  1900
# HR[5]        0.712   0.085    0.563    0.652    0.707    0.761    0.893 1.006   590
# beta[1]      5.491   0.223    5.262    5.381    5.463    5.552    5.892 1.075   590
# beta[2]      0.288   0.120    0.058    0.207    0.289    0.366    0.523 1.004   750
# beta[3]      0.354   0.102    0.157    0.285    0.352    0.424    0.554 1.002  4000
# beta[4]      0.301   0.100    0.109    0.235    0.300    0.368    0.497 1.002  2000
# beta[5]      0.229   0.079    0.078    0.180    0.230    0.281    0.379 1.009   440
# p.crit[2]    0.991   0.097    1.000    1.000    1.000    1.000    1.000 1.031  1600
# p.crit[3]    0.999   0.027    1.000    1.000    1.000    1.000    1.000 1.292  1300
# p.crit[4]    0.997   0.059    1.000    1.000    1.000    1.000    1.000 1.161   690
# p.crit[5]    0.995   0.070    1.000    1.000    1.000    1.000    1.000 1.160   490
# shape        1.510   0.084    1.345    1.474    1.517    1.557    1.637 1.054   510
# deviance  7294.526  49.407 7206.710 7261.987 7292.479 7324.320 7392.917 1.009  4000
# 
# For each parameter, n.eff is a crude measure of effective sample size,
# and Rhat is the potential scale reduction factor (at convergence, Rhat=1).
# 
# DIC info (using the rule, pD = var(deviance)/2)
# pD = 1221.0 and DIC = 8515.5
# DIC is an estimate of expected predictive error (lower deviance is better).

# Iterations = 1:97389
# Thinning interval = 97 
# Number of chains = 4 
# Sample size per chain = 1005 
# 
# 1. Empirical mean and standard deviation for each variable,
# plus standard error of the mean:
#   
#                Mean       SD  Naive SE Time-series SE
# AFT[2]       0.7549  0.09116 0.0014377      0.0030111
# AFT[3]       0.7058  0.07248 0.0011431      0.0026088
# AFT[4]       0.7436  0.07506 0.0011839      0.0026173
# AFT[5]       0.7976  0.06398 0.0010091      0.0026571
# HR[2]        0.6572  0.12096 0.0019079      0.0042032
# HR[3]        0.5931  0.09320 0.0014700      0.0037101
# HR[4]        0.6414  0.09794 0.0015447      0.0035681
# HR[5]        0.7115  0.08540 0.0013469      0.0038084
# beta[1]      5.4908  0.22345 0.0035242      0.0208103
# beta[2]      0.2884  0.12042 0.0018993      0.0039155
# beta[3]      0.3536  0.10225 0.0016128      0.0036670
# beta[4]      0.3013  0.10028 0.0015816      0.0034103
# beta[5]      0.2293  0.07872 0.0012416      0.0031069
# deviance  7294.5258 49.40739 0.7792538      2.0589413
# p.crit[2]    0.9905  0.09678 0.0015264      0.0035932
# p.crit[3]    0.9993  0.02731 0.0004308      0.0007651
# p.crit[4]    0.9965  0.05892 0.0009293      0.0026315
# p.crit[5]    0.9950  0.07037 0.0011098      0.0037671
# shape        1.5103  0.08353 0.0013174      0.0082829
# 
# 2. Quantiles for each variable:
#   
#                2.5%       25%       50%       75%     97.5%
# AFT[2]    5.927e-01    0.6934    0.7492    0.8130    0.9437
# AFT[3]    5.748e-01    0.6545    0.7031    0.7523    0.8548
# AFT[4]    6.086e-01    0.6918    0.7407    0.7904    0.8965
# AFT[5]    6.843e-01    0.7547    0.7949    0.8354    0.9250
# HR[2]     4.504e-01    0.5749    0.6464    0.7305    0.9190
# HR[3]     4.336e-01    0.5274    0.5874    0.6505    0.7970
# HR[4]     4.720e-01    0.5732    0.6356    0.7021    0.8478
# HR[5]     5.632e-01    0.6522    0.7068    0.7608    0.8930
# beta[1]   5.262e+00    5.3814    5.4626    5.5524    5.8919
# beta[2]   5.790e-02    0.2070    0.2887    0.3662    0.5231
# beta[3]   1.569e-01    0.2847    0.3523    0.4240    0.5537
# beta[4]   1.093e-01    0.2352    0.3001    0.3684    0.4966
# beta[5]   7.796e-02    0.1798    0.2296    0.2814    0.3794
# deviance  7.207e+03 7261.9869 7292.4791 7324.3202 7392.9169
# p.crit[2] 1.000e+00    1.0000    1.0000    1.0000    1.0000
# p.crit[3] 1.000e+00    1.0000    1.0000    1.0000    1.0000
# p.crit[4] 1.000e+00    1.0000    1.0000    1.0000    1.0000
# p.crit[5] 1.000e+00    1.0000    1.0000    1.0000    1.0000
# shape     1.345e+00    1.4742    1.5173    1.5571    1.6374

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
end.time - start.time # Time difference of 16.76311 hours


M1menIscheStroke_20200702 <- jagsfit
print(M1menIscheStroke_20200702)
summary(mcmcplots::as.mcmc.rjags(M1menIscheStroke_20200702))
save(M0menIscheStroke_20200702, M1menIscheStroke_20200702, file = "IscheStrokeMen.RData")


# Inference for Bugs model at "AgeIsche.weibull.model", fit using jags,
# 4 chains, each with 1e+05 iterations (first 2500 discarded), n.thin = 97
# n.sims = 4020 iterations saved
#            mu.vect sd.vect     2.5%      25%      50%      75%    97.5%  Rhat n.eff
# AFT[2]       0.835   0.085    0.678    0.778    0.832    0.889    1.010 1.004   640
# AFT[3]       0.792   0.066    0.671    0.746    0.790    0.836    0.930 1.009   300
# AFT[4]       0.821   0.069    0.694    0.773    0.819    0.863    0.964 1.005   630
# AFT[5]       0.740   0.046    0.655    0.710    0.738    0.768    0.835 1.006   440
# HR[2]        0.723   0.134    0.492    0.629    0.713    0.806    1.018 1.004   770
# HR[3]        0.655   0.102    0.479    0.582    0.649    0.721    0.879 1.010   280
# HR[4]        0.698   0.109    0.512    0.622    0.690    0.764    0.938 1.005   540
# HR[5]        0.576   0.067    0.459    0.531    0.570    0.615    0.723 1.006   490
# p.crit[2]    0.969   0.174    0.000    1.000    1.000    1.000    1.000 1.030   540
# p.crit[3]    0.998   0.047    1.000    1.000    1.000    1.000    1.000 1.022  4000
# p.crit[4]    0.990   0.100    1.000    1.000    1.000    1.000    1.000 1.023  2100
# p.crit[5]    1.000   0.000    1.000    1.000    1.000    1.000    1.000 1.000     1
# deviance  6403.727  52.849 6311.268 6371.579 6401.316 6433.091 6502.175 1.032   310
# 
# For each parameter, n.eff is a crude measure of effective sample size,
# and Rhat is the potential scale reduction factor (at convergence, Rhat=1).
# 
# DIC info (using the rule, pD = var(deviance)/2)
# pD = 1383.6 and DIC = 7787.4
# DIC is an estimate of expected predictive error (lower deviance is better).

# Iterations = 1:97389
# Thinning interval = 97 
# Number of chains = 4 
# Sample size per chain = 1005 
# 
# 1. Empirical mean and standard deviation for each variable,
# plus standard error of the mean:
#   
#                Mean       SD  Naive SE Time-series SE
# AFT[2]       0.8354  0.08503 0.0013411       0.002905
# AFT[3]       0.7923  0.06607 0.0010421       0.002182
# AFT[4]       0.8206  0.06887 0.0010863       0.002311
# AFT[5]       0.7397  0.04572 0.0007210       0.001804
# HR[2]        0.7232  0.13434 0.0021188       0.004551
# HR[3]        0.6550  0.10174 0.0016047       0.003541
# HR[4]        0.6983  0.10870 0.0017144       0.004067
# HR[5]        0.5762  0.06741 0.0010633       0.002864
# deviance  6403.7269 52.84895 0.8335340       2.602053
# p.crit[2]    0.9687  0.17427 0.0027485       0.004999
# p.crit[3]    0.9978  0.04727 0.0007455       0.001089
# p.crit[4]    0.9898  0.10049 0.0015849       0.003800
# p.crit[5]    1.0000  0.00000 0.0000000       0.000000
# 
# 2. Quantiles for each variable:
#   
#                2.5%       25%       50%       75%     97.5%
# AFT[2]       0.6781    0.7784    0.8320    0.8889    1.0100
# AFT[3]       0.6714    0.7459    0.7900    0.8364    0.9298
# AFT[4]       0.6940    0.7733    0.8186    0.8632    0.9643
# AFT[5]       0.6546    0.7103    0.7384    0.7681    0.8346
# HR[2]        0.4920    0.6293    0.7125    0.8057    1.0179
# HR[3]        0.4794    0.5824    0.6491    0.7205    0.8790
# HR[4]        0.5120    0.6223    0.6898    0.7639    0.9384
# HR[5]        0.4592    0.5311    0.5701    0.6155    0.7231
# deviance  6311.2677 6371.5789 6401.3162 6433.0914 6502.1750
# p.crit[2]    0.0000    1.0000    1.0000    1.0000    1.0000
# p.crit[3]    1.0000    1.0000    1.0000    1.0000    1.0000
# p.crit[4]    1.0000    1.0000    1.0000    1.0000    1.0000
# p.crit[5]    1.0000    1.0000    1.0000    1.0000    1.0000

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




M2Isch.weibull.model <- function() {
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
                         model.file=M2Isch.weibull.model)
end.time <- Sys.time()
end.time - start.time #Time difference of 2.743878 days


M2menIscheStroke_20200702 <- jagsfit
print(M2menIscheStroke_20200702)
summary(mcmcplots::as.mcmc.rjags(M2menIscheStroke_20200702))
save(M0menIscheStroke_20200702, 
     M1menIscheStroke_20200702, 
     M2menIscheStroke_20200702, file = "IscheStrokeMen.RData")

# Inference for Bugs model at "M2Isch.weibull.model", fit using jags,
# 4 chains, each with 1e+05 iterations (first 1500 discarded), n.thin = 98
# n.sims = 4020 iterations saved
#            mu.vect sd.vect     2.5%      25%      50%      75%    97.5%  Rhat n.eff
# AFT[2]       0.837   0.093    0.672    0.777    0.834    0.893    1.024 1.005   530
# AFT[3]       0.798   0.077    0.666    0.751    0.796    0.843    0.949 1.009   540
# AFT[4]       0.834   0.080    0.695    0.783    0.831    0.880    0.998 1.005   520
# AFT[5]       0.754   0.056    0.664    0.722    0.753    0.785    0.854 1.008   410
# HR[2]        0.731   0.140    0.496    0.631    0.719    0.817    1.040 1.004   680
# HR[3]        0.671   0.111    0.484    0.594    0.660    0.734    0.912 1.005   730
# HR[4]        0.725   0.120    0.523    0.640    0.714    0.794    0.996 1.004   800
# HR[5]        0.606   0.079    0.482    0.555    0.597    0.646    0.787 1.009   430
# p.crit[2]    0.961   0.194    0.000    1.000    1.000    1.000    1.000 1.013   990
# p.crit[3]    0.991   0.097    1.000    1.000    1.000    1.000    1.000 1.074   660
# p.crit[4]    0.975   0.155    1.000    1.000    1.000    1.000    1.000 1.015  1300
# p.crit[5]    0.998   0.045    1.000    1.000    1.000    1.000    1.000 1.295   500
# deviance  6397.166  83.808 6294.880 6354.341 6385.871 6418.905 6629.781 1.020   910
# 
# For each parameter, n.eff is a crude measure of effective sample size,
# and Rhat is the potential scale reduction factor (at convergence, Rhat=1).
# 
# DIC info (using the rule, pD = var(deviance)/2)
# pD = 3502.6 and DIC = 9899.8
# DIC is an estimate of expected predictive error (lower deviance is better).

#Iterations = 1:98393
# Thinning interval = 98 
# Number of chains = 4 
# Sample size per chain = 1005 
# 
# 1. Empirical mean and standard deviation for each variable,
# plus standard error of the mean:
#   
#                Mean       SD  Naive SE Time-series SE
# AFT[2]       0.8368  0.09312 0.0014687       0.003333
# AFT[3]       0.7979  0.07709 0.0012159       0.002909
# AFT[4]       0.8337  0.08000 0.0012618       0.003153
# AFT[5]       0.7543  0.05626 0.0008874       0.002461
# HR[2]        0.7314  0.13987 0.0022060       0.004901
# HR[3]        0.6707  0.11108 0.0017520       0.004738
# HR[4]        0.7245  0.12000 0.0018927       0.006061
# HR[5]        0.6058  0.07866 0.0012406       0.005193
# deviance  6397.1658 83.80783 1.3218177       8.655554
# p.crit[2]    0.9609  0.19375 0.0030558       0.008349
# p.crit[3]    0.9905  0.09678 0.0015264       0.004648
# p.crit[4]    0.9754  0.15500 0.0024447       0.008822
# p.crit[5]    0.9980  0.04457 0.0007030       0.001845
# 
# 2. Quantiles for each variable:
#   
#                2.5%       25%       50%       75%     97.5%
# AFT[2]       0.6724    0.7770    0.8344    0.8927    1.0238
# AFT[3]       0.6656    0.7513    0.7958    0.8427    0.9490
# AFT[4]       0.6948    0.7831    0.8307    0.8799    0.9978
# AFT[5]       0.6643    0.7224    0.7528    0.7853    0.8545
# HR[2]        0.4958    0.6310    0.7190    0.8169    1.0403
# HR[3]        0.4845    0.5939    0.6605    0.7341    0.9122
# HR[4]        0.5228    0.6401    0.7138    0.7944    0.9963
# HR[5]        0.4819    0.5547    0.5967    0.6456    0.7868
# deviance  6294.8797 6354.3411 6385.8713 6418.9045 6629.7809
# p.crit[2]    0.0000    1.0000    1.0000    1.0000    1.0000
# p.crit[3]    1.0000    1.0000    1.0000    1.0000    1.0000
# p.crit[4]    1.0000    1.0000    1.0000    1.0000    1.0000
# p.crit[5]    1.0000    1.0000    1.0000    1.0000    1.0000