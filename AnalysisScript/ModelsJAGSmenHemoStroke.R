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

# Inference for Bugs model at "M2.weibull.model", fit using jags,
# 3 chains, each with 1e+05 iterations (first 250 discarded), n.thin = 99
# n.sims = 3021 iterations saved
#             mu.vect sd.vect      2.5%       25%       50%       75%     97.5%  Rhat n.eff
# AFT[2]        1.006   0.124     0.847     0.942     0.995     1.053     1.203 1.009  1100
# AFT[3]        1.114   0.142     0.977     1.053     1.098     1.148     1.295 1.016   710
# AFT[4]        1.022   0.119     0.886     0.969     1.010     1.056     1.186 1.014  2900
# AFT[5]        0.972   0.097     0.877     0.934     0.963     0.994     1.104 1.015  1700
# HR[2]         1.010   0.168     0.752     0.903     0.992     1.092     1.357 1.004   660
# HR[3]         1.194   0.170     0.960     1.093     1.175     1.268     1.517 1.007   380
# HR[4]         1.033   0.150     0.813     0.947     1.018     1.097     1.314 1.004  1400
# HR[5]         0.950   0.117     0.798     0.889     0.938     0.990     1.175 1.006   610
# p.crit[2]     0.528   0.499     0.000     0.000     1.000     1.000     1.000 1.003   690
# p.crit[3]     0.064   0.245     0.000     0.000     0.000     0.000     1.000 1.007  1100
# p.crit[4]     0.444   0.497     0.000     0.000     0.000     1.000     1.000 1.003   920
# p.crit[5]     0.781   0.414     0.000     1.000     1.000     1.000     1.000 1.002  1000
# deviance  11949.549 126.233 11809.943 11893.025 11934.989 11981.600 12112.061 1.005  1800
# 
# For each parameter, n.eff is a crude measure of effective sample size,
# and Rhat is the potential scale reduction factor (at convergence, Rhat=1).
# 
# DIC info (using the rule, pD = var(deviance)/2)
# pD = 7964.3 and DIC = 19913.8
# DIC is an estimate of expected predictive error (lower deviance is better).




#######################Model 1 age-adj######################
#                                                          #
#                      Model 1 age-adj                     #
#                                                          #
#**********************************************************#
