# https://www.bioconductor.org/help/course-materials/2016/BioC2016/ConcurrentWorkshops4/Buros/weibull-survival-model.html


library(biostan)
stan_file <- system.file('stan', 'weibull_survival_null_model.stan', package =  'biostan')
biostan::print_stan_file(stan_file)
if (interactive())
  file.edit(stan_file)
print_stan_file(stan_file, section = 'data')
print_stan_file(stan_file, section = 'model')
print_stan_file(stan_file, section = 'transformed parameters')



sim_data <- function(alpha, mu, Nobs, Ncen) {
  observed_data <- data.frame(os_status = rep_len('DECEASED', Nobs),
                              os_months = rweibull(n = Nobs, alpha, exp(-(mu)/alpha)),
                              stringsAsFactors = F
  )
  
  censored_data <- data.frame(os_status = rep_len('LIVING', Ncen),
                              os_months = runif(Ncen) * rweibull(Ncen, alpha, exp(-(mu)/alpha)),
                              stringsAsFactors = F
  )
  
  return(observed_data %>% bind_rows(censored_data))
}

test_alpha <- 0.8
test_mu <- -3

## sample sizes from TCGA blca data
test_nobs <- 179 
test_ncen <- 230

## test these inputs for arbitrary values of alpha & mu
simulated_data <- 
  sim_data(alpha = test_alpha,
           mu = test_mu,
           Nobs = test_nobs,
           Ncen = test_ncen
  ) 
head(simulated_data)


print_stan_file(stan_file, section = 'data')


observed_data <- simulated_data %>%
  dplyr::filter(os_status == 'DECEASED')

censored_data <- simulated_data %>%
  dplyr::filter(os_status != 'DECEASED')

stan_data <- list(
  Nobs = nrow(observed_data),
  Ncen = nrow(censored_data),
  yobs = observed_data$os_months,
  ycen = censored_data$os_months
)
rm(censored_data)
rm(observed_data)
str(stan_data)


recover_simulated <- 
  rstan::stan(stan_file,
              data = stan_data,
              chains = 4,
              iter = 1000,
              seed = 1328025050
  )
print(recover_simulated)
