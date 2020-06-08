load("data/JACCmilkstroke.Rdata")


##%######################################################%##
#                                                          #
####          Example use kidney data in brms           ####
#                                                          #
##%######################################################%##


library(brms)
data("kidney")

head(kidney)
table(kidney$censored)
## performing surivival analysis using the "weibull" family
fit1 <- brm(time | cens(censored) ~ age + sex + disease, 
            data = kidney, family = weibull, inits = "0")
summary(fit1) 
plot(fit1)

##%######################################################%##
#                                                          #
####                 try use men data                   ####
#                                                          #
##%######################################################%##
library(tidyverse)
Mlkfre <- MData_men$Mlkfre

# is.censored <- is.na(t)
is.censored <- MData_men$Tot_Stroke != "I60_9"

# t <- larynx.dat$time
t <- if_else(!is.censored, MData_men$followpy, 0)
t <- na_if(t, 0)

Age <- MData_men$Age

c <- if_else(is.censored, MData_men$followpy, 0)


c[!is.censored] <- t[!is.censored] + 1

MILKdataMEN_brm <- data.frame(t = t,  c = c, Mlkfre = Mlkfre, is.censored = is.censored)


M0men20200608 <- brm(c | cens(is.censored) ~ Mlkfre, 
                     data = MILKdataMEN_brm, 
                     family = weibull, inits = "0")
