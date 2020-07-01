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
end.time - start.time

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
