library(biostan)
library(rstan)
library(bayesplot)
library(tidybayes)  
library(survival)
library(tidyverse)
library(survminer)

su_obj_men <- Surv(MData_men$followpy, MData_men$Tot_Stroke == "I60_9")

km_fit <- survfit(su_obj_men ~ Mlkfre, data = MData_men)
km_fit
ggsurvplot(km_fit,
           conf.int = TRUE,
           break.time.by = 4,
           risk.table = TRUE, 
           pval = TRUE,
           ylim = c(0.93, 1))

stan_weibull_survival_model_file <- system.file('stan', 'weibull_survival_model.stan', package =  'biostan')
biostan::print_stan_file(stan_weibull_survival_model_file)


stan_weibull_survival_model_data <-
  list(
    ## Number of event individuals
    Nobs = sum(MData_men$Tot_Stroke == "I60_9"),
    ## Number of censored individuals
    Ncen = sum(MData_men$Tot_Stroke != "I60_9"),
    ## Number of covariates
    M_bg = 1,
    ## Times for event individuals
    yobs = MData_men$followpy[MData_men$Tot_Stroke == "I60_9"],
    ## Times for censored individuals
    ycen = MData_men$followpy[MData_men$Tot_Stroke != "I60_9"],
    ## Covariates for event individuals as a matrix
    Xobs_bg = matrix(as.numeric(MData_men$MlkLogi == "Drinker")[MData_men$Tot_Stroke == "I60_9"]),
    ## Covariates for censored individuals as a matrix
    Xcen_bg = matrix(as.numeric(MData_men$MlkLogi == "Never")[MData_men$Tot_Stroke != "I60_9"])
  )
stan_weibull_survival_model_data

stan_weibull_survival_model_fit <-
  rstan::stan(file = stan_weibull_survival_model_file,
              data = stan_weibull_survival_model_data)
