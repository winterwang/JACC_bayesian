# stage <- larynx.dat$stage
Mlkfre <- MData_men$Mlkfre
MlkLogi <- as.numeric(MData_men$MlkLogi == "Drinker")

# is.censored <- is.na(t)
is.censored <- MData_men$Tot_Stroke != "I60_9"

# t <- larynx.dat$time
t <- if_else(!is.censored, MData_men$followpy, 0)
t <- na_if(t, 0)
t

# age <- larynx.dat$age
Agegrp <- MData_men$Agegrp
Age <- MData_men$Age


# c <- larynx.dat$cens_time
c <- if_else(is.censored, MData_men$followpy, 0)


c[!is.censored] <- t[!is.censored] + 1
# c
# t
# is.censored


larynx.weibull.model <- function() {
  for(i in 1:39386){
    sAge[i] <- (Age[i] - mean(Age[])) / sd(Age[])
    # sYr[i] <- (yr[i] - mean(yr[])) / sd(yr[])
    is.censored[i] ~ dinterval(t[i], c[i])
    t[i] ~ dweib(shape, lambda[i])
    lambda[i] <- exp(-mu[i] * shape)
    mu[i] <- beta[1] + beta[2]*equals(MlkLogi[i], 1) +
      # beta[3]*equals(Age[i], 3) + beta[4] * equals(stage[i], 4) +
      beta[3] *sAge[i] 
    
    ### calculate log-likelihoods
    y[i] <- ifelse(is.censored[i], c[i], t[i])
    loglik[i] <- log(ifelse(is.censored[i],
                            exp(-lambda[i] * (y[i] ^ shape)),
                            shape * lambda[i] * (y[i] ^ (shape - 1)) * exp(-lambda[i] * (y[i] ^ shape))))
  }
  
  ## priors for betas
  for(j in 1:3){
    beta[j] ~ dnorm(0, 0.001)
  }
  
  ### prior for shape
  shape ~ dgamma(.001, .001)
  
  ### Generated values
  AFT[2] <- exp(beta[2])
  HR[2] <- exp(shape * beta[2])
  p.crit[2] <- step(1 - HR[2])  
  AFT[3] <- exp(beta[3])
  HR[3] <- exp(shape * beta[3])
  p.crit[3] <- step(1 - HR[3])
}

# larynx.data <- c("t", "age", "yr", "c", "stage", "is.censored")
MILKdataMEN <- c("t", "Age", "c", "MlkLogi", "is.censored")
larynx.params <- c("beta", "AFT", "HR", "p.crit", "shape")
library(R2jags)

JACC.weibull.fit <- jags(data = MILKdataMEN,
                           parameters.to.save = larynx.params,
                           n.chains = 4,
                           n.iter = 200000,
                           n.burnin = 10000,
                           n.thin = 1,
                           model.file = larynx.weibull.model)
print(JACC.weibull.fit)
mcmcplots::traplot(JACC.weibull.fit)
Simulated <- coda::as.mcmc(JACC.weibull.fit)
# postsamples <- buildMCMC("*")
gelman.diag(Simulated)
gelman.plot(Simulated)

library(ggmcmc)
ggSample <- ggs(Simulated)
ggSample %>% 
  filter(Iteration >= 1001 & Parameter %in% c("beta[1]", "AFT[1]")) %>% 
  ggs_autocorrelation()

HR2 <- ggSample %>% 
  filter(Parameter == "HR[2]")

plot(density(HR2$value), main = "HR sample 10000", 
     ylab = "P(theta)", xlab = "theta", col = "red")
hist(Y$value, main = "y sample 10000", ylab = "P(Y)", 
     xlab = "y", col = "red", prob = TRUE)
