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
is.censored <- MData_men$HemoStroke == "I60_2"

# define followup time for events
t <- if_else(!is.censored, MData_men$followpy, 0)
t <- na_if(t, 0)
t


# define followup time for censored 
c <- if_else(is.censored, MData_men$followpy, 0)

c[!is.censored] <- t[!is.censored] + 1




jacc.weibull.model0 <- function() {
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
jagsfit <- jags.parallel(data=MILKdataMEN,  parameters.to.save = larynx.params, 
                         n.iter=100000, n.burnin=(5000/2), n.chains = 3, #n.thin = 100,
                         model.file=jacc.weibull.model0)
end.time <- Sys.time()
end.time - start.time 
M0menHemoStroke_20200618 <- jagsfit
print(M0menHemoStroke_20200618)
save.image(file = "data/JACCmilkstrokewithHemo.Rdata")
# save.image(file = "data/JACCmilkstroke.Rdata")




# Inference for Bugs model at "M2.weibull.model", fit using jags,
# 3 chains, each with 1e+05 iterations (first 250 discarded), n.thin = 99
# n.sims = 3021 iterations saved
# mu.vect sd.vect      2.5%       25%       50%       75%     97.5%  Rhat n.eff
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
