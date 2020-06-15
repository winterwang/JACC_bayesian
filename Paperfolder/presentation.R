load("data/JACCmilkstroke.Rdata")

summary(coda::as.mcmc(M0men_20200520))
summary(coda::as.mcmc(M1_agemen_20200520))
summary(coda::as.mcmc(M2men_20200521))
