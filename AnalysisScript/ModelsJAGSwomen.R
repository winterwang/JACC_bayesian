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
