library(tidyverse)
library(rstan)
library(tidybayes)
library(bayestestR)

## Load in the data ----
Gray_mrt <- read_delim("Gray_mrt.csv", ";", escape_double = FALSE,
                       col_types = cols(date_m = col_date(format = "%d.%m.%y"), date_r = col_date(format = "%d.%m.%y")),
                       trim_ws = TRUE)

# adding time difference (dt) and length increments (dL) to data
Gray_mrt <- Gray_mrt %>% 
  mutate(dt = as.numeric((difftime(date_r, date_m, units = "days") / 365)),
         dL = Lr - Lm)

priors <- read_csv2("Priors.csv")

## Bayesian Fabens model in stan ----

# Generating mean and sd for the prior probability distribution of Linf and K
mean_Linf <- mean(log(priors$Linf))
sd_Linf <- sd(log(priors$Linf))
mean_K <- mean(log(priors$K))
sd_K <- sd(log(priors$K))

# Putting data in a list for shipment to stan
data_list = list(
  Lm = Gray_mrt$Lm,
  Lr = Gray_mrt$Lr,
  dt = Gray_mrt$dt,
  N = length(Gray_mrt$Lm),
  mLinf = mean_Linf,
  sdLinf = sd_Linf,
  mK = mean_K,
  sdK = sd_K)

#Initial values for the MCMC chains
init_list <- lapply(1:4, function(i) {
  list(logLinf = log(rnorm(1, max(data_list$Lr), 1)), logK = log(runif(1,0.01,2)), sigma = runif(1,1,10))
})

# Compile and fit model to stan
stanfit_Fabens <- stan(file = "Bayesian_Fabens.stan",
                data = data_list,
                init = init_list,
                chains = 4,
                iter = 50000,
                warmup = 25000,
                thin = 10)

# Results and diagnostics for Fabens growth model ----
stanfit_Fabens

# View visual diagnostics for the Fabens model
traceplot(stanfit_Fabens, pars = c("Linf", "K", "sigma")) +
  ggtitle("Traceplot of the Fabens growth model")

stan_dens(stanfit_Fabens, pars = c("Linf", "K"))
