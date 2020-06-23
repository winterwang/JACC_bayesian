# jaccbayes


# stage <- larynx.dat$stage
Mlkfre <- as.numeric(MData_men$Mlkfre)
table(Mlkfre)
table(MData_men$Mlkfre) #1 = Never ; 2 = Mon1_2 ; 3 = Wek1_2 ; 4 = Wek3_4 ; 5 = Daily
MlkLogi <- as.numeric(MData_men$MlkLogi == "Drinker")

# is.censored <- is.na(t)
is.censored <- MData_men$Tot_Stroke != "I60_9"

# t <- larynx.dat$time
t <- if_else(!is.censored, MData_men$followpy, 0)
t <- na_if(t, 0)
t

# age <- larynx.dat$age
Agegrp <- as.numeric(MData_men$Agegrp)
table(Agegrp)
table(MData_men$Agegrp)  # 1 = [30,45) ; 2 = [45,55) ; 3 = [55,65) ; 4 =  [65,75) ; 5 = [75,80) 
Age <- MData_men$Age


# c <- larynx.dat$cens_time
c <- if_else(is.censored, MData_men$followpy, 0)


c[!is.censored] <- t[!is.censored] + 1
# c
# t
# is.censored


larynx.weibull.model <- function() {
  for(i in 1:39386){
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

# larynx.data <- c("t", "age", "yr", "c", "stage", "is.censored")
MILKdataMEN <- c("t",  "c", "Mlkfre", "is.censored")
larynx.params <- c("beta", "AFT", "HR", "p.crit", "shape")
library(R2jags)

# JACC.weibull.fit <- jags(data = MILKdataMEN,
#                          parameters.to.save = larynx.params,
#                          n.chains = 2,
#                          n.iter = 1000,
#                          n.burnin = 100,
#                          n.thin = 1,
#                          model.file = larynx.weibull.model)
start.time <- Sys.time()
jagsfit <- jags.parallel(data=MILKdataMEN,  parameters.to.save = larynx.params, 
                         n.iter=100000, n.burnin=(5000/2), n.chains = 3,
                         model.file=larynx.weibull.model)
end.time <- Sys.time()
end.time - start.time #Time difference of 18.74414 hours
M0men_20200520 <- jagsfit
# print(JACC.weibull.fit)
print(M0men_20200520)
save.image(file = "data/JACCmilkstroke.Rdata")

# 
# Iterations = 1:97389
# Thinning interval = 97 
# Number of chains = 3 
# Sample size per chain = 1005 
# 
# 1. Empirical mean and standard deviation for each variable,
# plus standard error of the mean:
#   
#                Mean       SD  Naive SE Time-series SE
# AFT[2]    9.276e-01  0.06579 0.0011982      0.0015899
# AFT[3]    8.345e-01  0.05350 0.0009744      0.0014405
# AFT[4]    8.491e-01  0.05445 0.0009917      0.0015845
# AFT[5]    9.312e-01  0.04423 0.0008056      0.0012693
# HR[2]     8.981e-01  0.09250 0.0016846      0.0022461
# HR[3]     7.701e-01  0.07135 0.0012995      0.0019555
# HR[4]     7.898e-01  0.07320 0.0013331      0.0021559
# HR[5]     9.024e-01  0.06216 0.0011320      0.0017923
# beta[1]   5.050e+00  0.06238 0.0011361      0.0034163
# beta[2]   7.768e-02  0.07096 0.0012923      0.0017130
# beta[3]   1.830e-01  0.06429 0.0011708      0.0016984
# beta[4]   1.656e-01  0.06412 0.0011678      0.0018635
# beta[5]   7.237e-02  0.04749 0.0008649      0.0013610
# deviance  1.649e+04 72.89115 1.3274896      1.5329732
# p.crit[2] 8.650e-01  0.34177 0.0062243      0.0068845
# p.crit[3] 9.990e-01  0.03153 0.0005743      0.0005743
# p.crit[4] 9.970e-01  0.05456 0.0009937      0.0011170
# p.crit[5] 9.347e-01  0.24717 0.0045014      0.0058620
# shape     1.453e+00  0.03510 0.0006393      0.0017323
# 
# 2. Quantiles for each variable:
#   
#                2.5%       25%       50%       75%     97.5%
# AFT[2]     8.055e-01 8.817e-01 9.264e-01 9.706e-01 1.060e+00
# AFT[3]     7.296e-01 7.983e-01 8.330e-01 8.710e-01 9.394e-01
# AFT[4]     7.443e-01 8.103e-01 8.464e-01 8.874e-01 9.596e-01
# AFT[5]     8.471e-01 9.007e-01 9.308e-01 9.605e-01 1.020e+00
# HR[2]      7.295e-01 8.336e-01 8.951e-01 9.577e-01 1.088e+00
# HR[3]      6.337e-01 7.210e-01 7.670e-01 8.180e-01 9.127e-01
# HR[4]      6.551e-01 7.375e-01 7.848e-01 8.410e-01 9.426e-01
# HR[5]      7.863e-01 8.587e-01 9.006e-01 9.429e-01 1.029e+00
# beta[1]    4.931e+00 5.007e+00 5.048e+00 5.092e+00 5.173e+00
# beta[2]   -5.849e-02 2.989e-02 7.640e-02 1.259e-01 2.163e-01
# beta[3]    6.250e-02 1.381e-01 1.827e-01 2.253e-01 3.152e-01
# beta[4]    4.121e-02 1.195e-01 1.668e-01 2.103e-01 2.953e-01
# beta[5]   -1.986e-02 4.032e-02 7.174e-02 1.045e-01 1.659e-01
# deviance   1.635e+04 1.644e+04 1.649e+04 1.654e+04 1.663e+04
# p.crit[2]  0.000e+00 1.000e+00 1.000e+00 1.000e+00 1.000e+00
# p.crit[3]  1.000e+00 1.000e+00 1.000e+00 1.000e+00 1.000e+00
# p.crit[4]  1.000e+00 1.000e+00 1.000e+00 1.000e+00 1.000e+00
# p.crit[5]  0.000e+00 1.000e+00 1.000e+00 1.000e+00 1.000e+00
# shape      1.385e+00 1.428e+00 1.452e+00 1.476e+00 1.523e+00



# Inference for Bugs model at "larynx.weibull.model", fit using jags,
# 3 chains, each with 1e+05 iterations (first 2500 discarded), n.thin = 97
# n.sims = 3015 iterations saved
#             mu.vect sd.vect      2.5%       25%       50%       75%     97.5%  Rhat n.eff
# AFT[2]        0.928   0.066     0.806     0.882     0.926     0.971     1.060 1.003   960
# AFT[3]        0.834   0.054     0.730     0.798     0.833     0.871     0.939 1.004   610
# AFT[4]        0.849   0.054     0.744     0.810     0.846     0.887     0.960 1.004   630
# AFT[5]        0.931   0.044     0.847     0.901     0.931     0.960     1.020 1.004   510
# HR[2]         0.898   0.093     0.730     0.834     0.895     0.958     1.088 1.003   960
# HR[3]         0.770   0.071     0.634     0.721     0.767     0.818     0.913 1.004   630
# HR[4]         0.790   0.073     0.655     0.738     0.785     0.841     0.943 1.004   640
# HR[5]         0.902   0.062     0.786     0.859     0.901     0.943     1.029 1.004   510
# beta[1]       5.050   0.062     4.931     5.007     5.048     5.092     5.173 1.002  1500
# beta[2]       0.078   0.071    -0.058     0.030     0.076     0.126     0.216 1.003   960
# beta[3]       0.183   0.064     0.062     0.138     0.183     0.225     0.315 1.004   610
# beta[4]       0.166   0.064     0.041     0.119     0.167     0.210     0.295 1.004   630
# beta[5]       0.072   0.047    -0.020     0.040     0.072     0.105     0.166 1.004   510
# p.crit[2]     0.865   0.342     0.000     1.000     1.000     1.000     1.000 1.004  1000
# p.crit[3]     0.999   0.032     1.000     1.000     1.000     1.000     1.000 1.134  3000
# p.crit[4]     0.997   0.055     1.000     1.000     1.000     1.000     1.000 1.067  2300
# p.crit[5]     0.935   0.247     0.000     1.000     1.000     1.000     1.000 1.002  3000
# shape         1.453   0.035     1.385     1.428     1.452     1.476     1.523 1.001  3000
# deviance  16489.564  72.891 16353.174 16439.538 16488.632 16540.869 16634.684 1.001  3000
# 
# For each parameter, n.eff is a crude measure of effective sample size,
# and Rhat is the potential scale reduction factor (at convergence, Rhat=1).
# 
# DIC info (using the rule, pD = var(deviance)/2)
# pD = 2658.2 and DIC = 19147.8
# DIC is an estimate of expected predictive error (lower deviance is better).

# mcmcplots::traplot(JACC.weibull.fit, c("beta[5]", "HR[5]"))
mcmcplots::traplot(jagsfit, c("beta[5]", "HR[5]"))

# samplesHistory("*", mfrow = c(3,1), beg = 501, ask = FALSE)
Simulated <- coda::as.mcmc(JACC.weibull.fit)
library(ggmcmc)

ggSample <- ggs(Simulated)
ggSample %>% 
  filter(Iteration >= 500 & Parameter %in% c("beta[5]", "HR[5]")) %>% 
  ggs_traceplot()

# update(JACC.weibull.fit, n.iter = 1000)
# autojags(JACC.weibull.fit, n.iter = 5000, Rhat = 1.05, n.update = 2)

# postsamples <- buildMCMC("*")
# gelman.diag(Simulated)
# gelman.plot(Simulated)
Simulated <- coda::as.mcmc(JACC.weibull.fit)

ggSample <- ggs(Simulated)
ggSample %>% 
  filter(Iteration >= 100 & Parameter %in% c("beta[5]", "AFT[5]")) %>% 
  ggs_autocorrelation()

HR5 <- ggSample %>% 
  filter(Parameter == "HR[5]")

plot(density(HR5$value), main = "HR sample 10000", 
     ylab = "P(theta)", xlab = "theta", col = "red")
# hist(Y$value, main = "y sample 10000", ylab = "P(Y)", 
#      xlab = "y", col = "red", prob = TRUE)

##%######################################################%##
#                                                          #
####             age adjusted model in MEN              ####
#                                                          #
##%######################################################%##

Age.weibull.model <- function() {
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
                         n.iter=100000, n.burnin=(5000/2), n.chains = 3,
                         model.file=Age.weibull.model)
end.time <- Sys.time()
end.time - start.time #Time difference of 1.093635 days

# Iterations = 1:97389
# Thinning interval = 97 
# Number of chains = 3 
# Sample size per chain = 1005 
# 
# 1. Empirical mean and standard deviation for each variable,
# plus standard error of the mean:
#   
#                Mean       SD  Naive SE Time-series SE
# AFT[2]    9.862e-01  0.06272 0.0011422       0.001589
# AFT[3]    9.025e-01  0.04964 0.0009040       0.001393
# AFT[4]    9.122e-01  0.04965 0.0009042       0.001329
# AFT[5]    8.476e-01  0.03511 0.0006394       0.001054
# HR[2]     9.790e-01  0.10599 0.0019304       0.002686
# HR[3]     8.411e-01  0.07881 0.0014353       0.002236
# HR[4]     8.565e-01  0.07926 0.0014435       0.002148
# HR[5]     7.552e-01  0.05329 0.0009705       0.001635
# deviance  1.488e+04 71.25668 1.2977228       1.623196
# p.crit[2] 5.867e-01  0.49250 0.0089694       0.011160
# p.crit[3] 9.731e-01  0.16172 0.0029452       0.003207
# p.crit[4] 9.605e-01  0.19474 0.0035466       0.003689
# p.crit[5] 1.000e+00  0.00000 0.0000000       0.000000
# 
# 2. Quantiles for each variable:
#   
#                2.5%       25%       50%       75%     97.5%
# AFT[2]    8.686e-01 9.412e-01 9.855e-01 1.029e+00 1.109e+00
# AFT[3]    8.079e-01 8.687e-01 9.007e-01 9.357e-01 1.001e+00
# AFT[4]    8.181e-01 8.787e-01 9.117e-01 9.452e-01 1.011e+00
# AFT[5]    7.816e-01 8.244e-01 8.466e-01 8.702e-01 9.198e-01
# HR[2]     7.865e-01 9.017e-01 9.754e-01 1.051e+00 1.193e+00
# HR[3]     6.977e-01 7.864e-01 8.366e-01 8.932e-01 1.002e+00
# HR[4]     7.108e-01 8.020e-01 8.533e-01 9.085e-01 1.020e+00
# HR[5]     6.566e-01 7.192e-01 7.522e-01 7.891e-01 8.666e-01
# deviance  1.475e+04 1.483e+04 1.488e+04 1.493e+04 1.502e+04
# p.crit[2] 0.000e+00 0.000e+00 1.000e+00 1.000e+00 1.000e+00
# p.crit[3] 0.000e+00 1.000e+00 1.000e+00 1.000e+00 1.000e+00
# p.crit[4] 0.000e+00 1.000e+00 1.000e+00 1.000e+00 1.000e+00
# p.crit[5] 1.000e+00 1.000e+00 1.000e+00 1.000e+00 1.000e+00

# Inference for Bugs model at "Age.weibull.model", fit using jags,
# 3 chains, each with 1e+05 iterations (first 2500 discarded), n.thin = 97
# n.sims = 3015 iterations saved
#             mu.vect sd.vect      2.5%       25%       50%       75%     97.5%  Rhat n.eff
# AFT[2]        0.986   0.063     0.869     0.941     0.986     1.029     1.109 1.001  2700
# AFT[3]        0.902   0.050     0.808     0.869     0.901     0.936     1.001 1.002  1600
# AFT[4]        0.912   0.050     0.818     0.879     0.912     0.945     1.011 1.003   870
# AFT[5]        0.848   0.035     0.782     0.824     0.847     0.870     0.920 1.003   850
# HR[2]         0.979   0.106     0.787     0.902     0.975     1.051     1.193 1.001  2800
# HR[3]         0.841   0.079     0.698     0.786     0.837     0.893     1.002 1.002  1200
# HR[4]         0.856   0.079     0.711     0.802     0.853     0.909     1.020 1.004   670
# HR[5]         0.755   0.053     0.657     0.719     0.752     0.789     0.867 1.005   490
# p.crit[2]     0.587   0.493     0.000     0.000     1.000     1.000     1.000 1.001  3000
# p.crit[3]     0.973   0.162     0.000     1.000     1.000     1.000     1.000 1.004  3000
# p.crit[4]     0.961   0.195     0.000     1.000     1.000     1.000     1.000 1.001  3000
# p.crit[5]     1.000   0.000     1.000     1.000     1.000     1.000     1.000 1.000     1
# deviance  14883.855  71.257 14747.919 14834.265 14883.382 14931.578 15022.621 1.003   880
# 
# For each parameter, n.eff is a crude measure of effective sample size,
# and Rhat is the potential scale reduction factor (at convergence, Rhat=1).
# 
# DIC info (using the rule, pD = var(deviance)/2)
# pD = 2534.7 and DIC = 17418.5
# DIC is an estimate of expected predictive error (lower deviance is better).

M1_agemen_20200520 <- jagsfit
save.image(file = "data/JACCmilkstroke.Rdata")

print(M1_agemen_20200520)
# mcmcplots::traplot(JACC.weibull.fit, c("beta[5]", "HR[5]"))
mcmcplots::traplot(jagsfit, c("AFT[5]", "HR[5]"))

##%######################################################%##
#                                                          #
####       Model 3, smking, alcohol, bmi, DM, HT,       ####
####  KID, LIV, Excercise, Slep, Spi, Cofe, Educ, Gret  ####
#                                                          #
##%######################################################%##


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


M2.weibull.model <- function() {
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
                         n.iter=100000, n.burnin=(3000/2), n.chains = 3,
                         model.file=M2.weibull.model)
end.time <- Sys.time()
end.time - start.time # Time difference of 4.046775 days
M2men_20200521 <- jagsfit
print(M2men_20200521)
save.image(file = "data/JACCmilkstroke.Rdata")


# Iterations = 1:98393
# Thinning interval = 98 
# Number of chains = 3 
# Sample size per chain = 1005 
# 
# 1. Empirical mean and standard deviation for each variable,
# plus standard error of the mean:
#   
#               Mean       SD  Naive SE Time-series SE
# AFT[2]    1.003e+00  0.07146 0.0013014       0.002869
# AFT[3]    9.195e-01  0.05815 0.0010591       0.002354
# AFT[4]    9.366e-01  0.06082 0.0011076       0.002560
# AFT[5]    8.782e-01  0.04600 0.0008378       0.002045
# HR[2]     1.005e+00  0.11614 0.0021152       0.004074
# HR[3]     8.679e-01  0.09041 0.0016465       0.003600
# HR[4]     8.953e-01  0.09461 0.0017231       0.003810
# HR[5]     8.022e-01  0.06951 0.0012660       0.003083
# deviance  1.484e+04 77.87825 1.4183144       2.532632
# p.crit[2] 5.061e-01  0.50005 0.0091068       0.012294
# p.crit[3] 9.373e-01  0.24244 0.0044153       0.007809
# p.crit[4] 8.962e-01  0.30507 0.0055559       0.008950
# p.crit[5] 9.904e-01  0.09762 0.0017778       0.005731
# 
# 2. Quantiles for each variable:
#   
#                2.5%       25%       50%       75%     97.5%
# AFT[2]    8.830e-01 9.562e-01 9.987e-01 1.043e+00 1.136e+00
# AFT[3]    8.238e-01 8.820e-01 9.161e-01 9.498e-01 1.032e+00
# AFT[4]    8.379e-01 8.988e-01 9.319e-01 9.685e-01 1.047e+00
# AFT[5]    8.073e-01 8.516e-01 8.750e-01 9.004e-01 9.610e-01
# HR[2]     8.070e-01 9.267e-01 9.979e-01 1.075e+00 1.236e+00
# HR[3]     7.177e-01 8.078e-01 8.613e-01 9.158e-01 1.053e+00
# HR[4]     7.383e-01 8.336e-01 8.866e-01 9.466e-01 1.082e+00
# HR[5]     6.946e-01 7.599e-01 7.968e-01 8.355e-01 9.346e-01
# deviance  1.470e+04 1.479e+04 1.483e+04 1.488e+04 1.498e+04
# p.crit[2] 0.000e+00 0.000e+00 1.000e+00 1.000e+00 1.000e+00
# p.crit[3] 0.000e+00 1.000e+00 1.000e+00 1.000e+00 1.000e+00
# p.crit[4] 0.000e+00 1.000e+00 1.000e+00 1.000e+00 1.000e+00
# p.crit[5] 1.000e+00 1.000e+00 1.000e+00 1.000e+00 1.000e+00

# Inference for Bugs model at "M2.weibull.model", fit using jags,
# 3 chains, each with 1e+05 iterations (first 1500 discarded), n.thin = 98
# n.sims = 3015 iterations saved
#             mu.vect sd.vect      2.5%       25%       50%       75%     97.5%  Rhat n.eff
# AFT[2]        1.003   0.071     0.883     0.956     0.999     1.043     1.136 1.004   870
# AFT[3]        0.919   0.058     0.824     0.882     0.916     0.950     1.032 1.006   520
# AFT[4]        0.937   0.061     0.838     0.899     0.932     0.969     1.047 1.004   680
# AFT[5]        0.878   0.046     0.807     0.852     0.875     0.900     0.961 1.004  1100
# HR[2]         1.005   0.116     0.807     0.927     0.998     1.075     1.236 1.003   920
# HR[3]         0.868   0.090     0.718     0.808     0.861     0.916     1.053 1.005   520
# HR[4]         0.895   0.095     0.738     0.834     0.887     0.947     1.082 1.004   660
# HR[5]         0.802   0.070     0.695     0.760     0.797     0.836     0.935 1.002  1000
# p.crit[2]     0.506   0.500     0.000     0.000     1.000     1.000     1.000 1.001  3000
# p.crit[3]     0.937   0.242     0.000     1.000     1.000     1.000     1.000 1.009   900
# p.crit[4]     0.896   0.305     0.000     1.000     1.000     1.000     1.000 1.001  3000
# p.crit[5]     0.990   0.098     1.000     1.000     1.000     1.000     1.000 1.009  3000
# deviance  14836.944  77.878 14702.544 14787.079 14833.169 14881.968 14983.800 1.002  1800
# 
# For each parameter, n.eff is a crude measure of effective sample size,
# and Rhat is the potential scale reduction factor (at convergence, Rhat=1).
# 
# DIC info (using the rule, pD = var(deviance)/2)
# pD = 3031.2 and DIC = 17868.1
# DIC is an estimate of expected predictive error (lower deviance is better).
# =======
# jaccbayes


# stage <- larynx.dat$stage
Mlkfre <- as.numeric(MData_men$Mlkfre)
table(Mlkfre)
table(MData_men$Mlkfre) #1 = Never ; 2 = Mon1_2 ; 3 = Wek1_2 ; 4 = Wek3_4 ; 5 = Daily
MlkLogi <- as.numeric(MData_men$MlkLogi == "Drinker")

# is.censored <- is.na(t)
is.censored <- MData_men$Tot_Stroke != "I60_9"

# t <- larynx.dat$time
t <- if_else(!is.censored, MData_men$followpy, 0)
t <- na_if(t, 0)
t

# age <- larynx.dat$age
Agegrp <- as.numeric(MData_men$Agegrp)
table(Agegrp)
table(MData_men$Agegrp)  # 1 = [30,45) ; 2 = [45,55) ; 3 = [55,65) ; 4 =  [65,75) ; 5 = [75,80) 
Age <- MData_men$Age


# c <- larynx.dat$cens_time
c <- if_else(is.censored, MData_men$followpy, 0)


c[!is.censored] <- t[!is.censored] + 1
# c
# t
# is.censored


larynx.weibull.model <- function() {
  for(i in 1:39386){
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

# larynx.data <- c("t", "age", "yr", "c", "stage", "is.censored")
MILKdataMEN <- c("t",  "c", "Mlkfre", "is.censored")
larynx.params <- c("beta", "AFT", "HR", "p.crit", "shape")
library(R2jags)

# JACC.weibull.fit <- jags(data = MILKdataMEN,
#                          parameters.to.save = larynx.params,
#                          n.chains = 2,
#                          n.iter = 1000,
#                          n.burnin = 100,
#                          n.thin = 1,
#                          model.file = larynx.weibull.model)

jagsfit <- jags.parallel(data=MILKdataMEN,  parameters.to.save = larynx.params, 
                         n.iter=10000, n.burnin=(5000/2), n.chains = 3,
                         model.file=larynx.weibull.model)

# print(JACC.weibull.fit)
print(jagsfit)
# mcmcplots::traplot(JACC.weibull.fit, c("beta[5]", "HR[5]"))
mcmcplots::traplot(jagsfit, c("beta[5]", "HR[5]"))

# samplesHistory("*", mfrow = c(3,1), beg = 501, ask = FALSE)
Simulated <- coda::as.mcmc(JACC.weibull.fit)
library(ggmcmc)

ggSample <- ggs(Simulated)
ggSample %>% 
  filter(Iteration >= 500 & Parameter %in% c("beta[5]", "HR[5]")) %>% 
  ggs_traceplot()

# update(JACC.weibull.fit, n.iter = 1000)
# autojags(JACC.weibull.fit, n.iter = 5000, Rhat = 1.05, n.update = 2)

# postsamples <- buildMCMC("*")
# gelman.diag(Simulated)
# gelman.plot(Simulated)
Simulated <- coda::as.mcmc(JACC.weibull.fit)

ggSample <- ggs(Simulated)
ggSample %>% 
  filter(Iteration >= 100 & Parameter %in% c("beta[5]", "AFT[5]")) %>% 
  ggs_autocorrelation()

HR5 <- ggSample %>% 
  filter(Parameter == "HR[5]")

plot(density(HR5$value), main = "HR sample 10000", 
     ylab = "P(theta)", xlab = "theta", col = "red")
# hist(Y$value, main = "y sample 10000", ylab = "P(Y)", 
#      xlab = "y", col = "red", prob = TRUE)

##%######################################################%##
#                                                          #
####             age adjusted model in MEN              ####
#                                                          #
##%######################################################%##

Age.weibull.model <- function() {
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
jagsfit <- jags.parallel(data=MILKdataMEN,  parameters.to.save = age.params, 
                         n.iter=10000, n.burnin=(5000/2), n.chains = 3,
                         model.file=Age.weibull.model)
print(jagsfit)
# mcmcplots::traplot(JACC.weibull.fit, c("beta[5]", "HR[5]"))
mcmcplots::traplot(jagsfit, c("AFT[5]", "HR[5]"))

##%######################################################%##
#                                                          #
####       Model 3, smking, alcohol, bmi, DM, HT,       ####
####  KID, LIV, Excercise, Slep, Spi, Cofe, Educ, Gret  ####
#                                                          #
##%######################################################%##


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


M2.weibull.model <- function() {
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
      beta[16] * equals(BMIgrp[i], 4) + beta[17] * equals(BMIgrp[i], 5) + 
     beta[18] * equals(DM_hist[i], 2) + beta[19] * equals(DM_hist[i], 3) +
    beta[20] * equals(HT_hist[i], 2) + beta[21] * equals(HT_hist[i], 3) +
    beta[22] * equals(KID_hist[i], 2) + beta[23] * equals(KID_hist[i], 3)  +
    beta[24] * equals(LIV_hist[i], 2) + beta[25] * equals(LIV_hist[i], 3)  +
    beta[26] * equals(Exercise[i], 2) + beta[27] * equals(Exercise[i], 3)  +
    beta[28] * equals(Slepgrp[i], 2) + beta[29] * equals(Slepgrp[i], 3) +
    beta[30] * equals(Slepgrp[i], 4) + beta[31] * equals(Slepgrp[i], 5) +
    beta[32] * equals(Spi[i], 2) + beta[33] * equals(Spi[i], 3) +
    beta[34] * equals(Spi[i], 4) + beta[35] * equals(Spi[i], 5) +
    beta[36] * equals(Cofe[i], 2) + beta[37] * equals(Cofe[i], 3) +
    beta[38] * equals(Cofe[i], 4) + beta[39] * equals(Educgrp[i], 2) +
    beta[40] * equals(Educgrp[i], 3) +  beta[41] * equals(Gretea[i], 2) +
    beta[42] * equals(Gretea[i], 3) +  beta[43] * equals(Gretea[i], 4) +
    beta[44] * equals(Fru[i], 2) + beta[45] * equals(Fru[i], 3) +
    beta[46] * equals(Fru[i], 4) + beta[47] * equals(Fru[i], 5) #+
    # beta[45] * equals(BMIgrp[i], 3) + beta[46] * equals(BMIgrp[i], 4) + 
    # beta[47] * equals(BMIgrp[i], 5) 
  }
  
  ## priors for betas
  for(j in 1:47){
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
                 "Alc_Fre", "BMIgrp", 
                 "DM_hist", "HT_hist", "KID_hist",
                 "LIV_hist", "Exercise", "Slepgrp", "Spi", "Cofe",
                 "Educgrp", "Gretea", "Fru"
)

M2.params <- c("AFT", "HR", "p.crit")
M2menfit <- jags.parallel(data=MILKdataMEN,  parameters.to.save = M2.params, 
                         n.iter=10000, n.burnin=(5000/2), n.chains = 3,
                         model.file=M2.weibull.model)
print(M2menfit)
