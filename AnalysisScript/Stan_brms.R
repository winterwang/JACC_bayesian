############################################################
#                                                          #
#              Example from the brms package               #
#                                                          #
############################################################

library(brms)
data("kidney")
fit1 <- brm(time | cens(censored) ~ age + sex + disease, 
            data = kidney, family = weibull, inits = "0")
summary(fit1) 
plot(fit1)

# ## adding random intercepts over patients
# fit2 <- brm(time | cens(censored) ~ age + sex + disease + (1|patient), 
#             data = kidney, family = weibull(), inits = "0",
#             prior = set_prior("cauchy(0,2)", class = "sd"))
# summary(fit2) 
# plot(fit2)         

############################################################
#                                                          #
#              do the crude model using brms               #
#                                                          #
############################################################

load("data/JACCmilkstroke.Rdata")

# define exposure
Mlkfre <- MData_men$Mlkfre 
table(Mlkfre)

# define event
is.censored <- as.integer(MData_men$Tot_Stroke == "I60_9")

# define time of followup
t <- MData_men$followpy

CrudeData <- data.frame(t, is.censored, Mlkfre)
# fit a crude weibull model
M0 <- brm(t | cens(is.censored) ~ Mlkfre, 
            data = CrudeData, family = weibull, inits = "0",
          cores = parallel::detectCores()-1)
summary(M0) 
plot(M0)

stancode(M0)
M0data <- standata(M0)

library(rstan)
M0_0 <- stan(file = "AnalysisScript/M0StanMen.stan", 
             data = M0data, 
             seed = 1234, 
             pars = c("HR", "AFT", "b", "shape"), 
             options(mc.cores = 8), 
             chains = 4, 
             iter = 10000)


############################################################
#                                                          #
#               Fit the age adjusted model                 #
#                                                          #
############################################################


Age <- MData_men$Age

AgeData <- data.frame(t, is.censored, Mlkfre, Age)


M1 <- brm(t | cens(is.censored) ~ Mlkfre + Age, 
          data = AgeData, family = weibull, inits = "0",
          cores = parallel::detectCores()-1)
summary(M1) 
plot(M1)