# example from the internet  https://jfiksel.github.io/2017-05-24-waic_aft_models_jags/

library(data.table)
larynx.url <- "http://people.oregonstate.edu/~calverta/BIDA/Chapter13/Larynx-Cancer-Data.txt"
larynx.dat <- fread(larynx.url)
larynx.dat <- larynx.dat[1:90, c(1, 2, 4, 5, 6)]
colnames(larynx.dat) <- c("stage", "time", "age", "yr", "cens_time")


stage <- larynx.dat$stage
t <- larynx.dat$time
age <- larynx.dat$age
yr <- larynx.dat$yr
c <- larynx.dat$cens_time
is.censored <- is.na(t)
c[!is.censored] <- t[!is.censored] + 1
c
t
is.censored


larynx.weibull.model <- function() {
  for(i in 1:90){
    sAge[i] <- (age[i] - mean(age[])) / sd(age[])
    sYr[i] <- (yr[i] - mean(yr[])) / sd(yr[])
    is.censored[i] ~ dinterval(t[i], c[i])
    t[i] ~ dweib(shape, lambda[i])
    lambda[i] <- exp(-mu[i] * shape)
    mu[i] <- beta[1] + beta[2]*equals(stage[i], 2) +
      beta[3]*equals(stage[i], 3) + beta[4] * equals(stage[i], 4) +
      beta[5] *sAge[i] + beta[6] *sYr[i]
    
    
    ### calculate log-likelihoods
    y[i] <- ifelse(is.censored[i], c[i], t[i])
    loglik[i] <- log(ifelse(is.censored[i],
                            exp(-lambda[i] * (y[i] ^ shape)),
                            shape * lambda[i] * (y[i] ^ (shape - 1)) * exp(-lambda[i] * (y[i] ^ shape))))
  }
  
  ## priors for betas
  for(j in 1:6){
    beta[j] ~ dnorm(0, 0.001)
  }
  
  ### prior for shape
  shape ~ dgamma(.001, .001)
  
  ### Generated values
  AFT[2] <- exp(beta[2])
  HR[2] <- exp(shape * beta[2])
  p.crit[2] <- step(1 - HR[2])
}

larynx.data <- c("t", "age", "yr", "c", "stage", "is.censored")
larynx.params <- c("beta", "AFT", "HR", "p.crit", "shape")
library(R2jags)

larynx.weibull.fit <- jags(data = larynx.data,
                           parameters.to.save = larynx.params,
                           n.chains = 2,
                           n.iter = 10000,
                           n.burnin = 1000,
                           n.thin = 3,
                           model.file = larynx.weibull.model)
print(larynx.weibull.fit)
