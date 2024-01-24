library(tidyverse)
library(TropFishR)
library(rstan)
library(tidybayes)
library(bayesplot)

## Load in the data ----
# Age and Length data
AL_gray <- read_delim("AgeLen_com.csv", ";", escape_double = FALSE, locale = locale(decimal_mark = ","), 
                      trim_ws = TRUE)

# Mark and recapture data and adding time difference (dt) and Length increments (dL)
Gray_mrt <- read_delim("Gray_mrt.csv", ";", escape_double = FALSE,
                       col_types = cols(date_m = col_date(format = "%d.%m.%y"), date_r = col_date(format = "%d.%m.%y")),
                       trim_ws = TRUE) %>%
  mutate(dt = as.numeric((difftime(date_r, date_m, units = "days") / 365)),
         dL = Lr - Lm)

# Prior information
priors <- read_csv2("Priors.csv")

# Length frequency data for ELEFAN
lfq_data <- read_delim("lfq_dat.csv", ";", escape_double = FALSE, col_types = cols(date = col_date(format = "%d.%m.%y"), 
                       date_season = col_date(format = "%d.%m.%y")), trim_ws = TRUE)

# Generating mean and sd for the prior probability distribution of Linf, K, and t0
# -> shared for typical vB and Fabens
mean_Linf <- mean(log(priors$Linf))
sd_Linf <- sd(log(priors$Linf))
mean_K <- mean(log(priors$K))
sd_K <- sd(log(priors$K))
mean_t0 <- mean(priors$t0)
sd_t0 <- sd(priors$t0)

## Bayesian von Bertalanffy growth model ----
# Data list for shipment to stan
data_listVB <- list(
  Age = AL_gray$age,
  TL = AL_gray$tl,
  N = length(AL_gray$age),
  mLinf = mean_Linf,
  sdLinf = sd_Linf,
  mK = mean_K,
  sdK = sd_K,
  mt0 = mean_t0,
  sdt0 = sd_t0)
# Initial values for MCMC chains
init_listVB <- lapply(1:4, function(i) {
  list(logLinf = log(rnorm(1, max(data_listVB$TL), 1)),
       logK = log(runif(1,0.01,5)),
       t0 = rnorm(1,1,1),
       sigma = runif(1,1,10))
})

# Compile nd fit model to stan
satnfit_vBert <- stan(file = "Bayesian_vBGM.stan",
                      data = data_listVB,
                      init = init_listVB,
                      chains = 4,
                      iter = 50000,
                      warmup = 25000,
                      thin = 10)

# Results and diagnostics for von Bertalanffy growth model ----
stanfit_vBert

# View visual diagnostics for the Fabens model
traceplot(stanfit_Fabens, pars = c("Linf", "K", "t0", "sigma")) +
  ggtitle("Traceplot of the von Bertalanffy growth model")

stan_dens(stanfit_Fabens, pars = c("Linf", "K", "t0", "sigma"))

# Get Credible Intervals for Linf and K
post_vB <- as.data.frame(stanfit_vBert) %>% 
  summarise(Lci_lo = quantile(Linf, probs = .05),
            Lci_hi = quantile(Linf, probs = .95),
            Kci_lo = quantile(K, probs = .05),
            Kci_hi = quantile(K, probs = .95),
            t0ci_lo = quantile(t0, probs = .05),
            t0ci_hi = quantile(t0, probs = .95))

## Bayesian Fabens growth model in stan ----

# Putting data in a list for shipment to stan
data_listF = list(
  Lm = Gray_mrt$Lm,
  Lr = Gray_mrt$Lr,
  dt = Gray_mrt$dt,
  N = length(Gray_mrt$Lm),
  mLinf = mean_Linf,
  sdLinf = sd_Linf,
  mK = mean_K,
  sdK = sd_K)

#Initial values for the MCMC chains
init_listF <- lapply(1:4, function(i) {
  list(logLinf = log(rnorm(1, max(data_list$Lr), 1)),
       logK = log(runif(1, mean(log(data_list$Lr) - log(data_list$Lm)) / data_list$dt, 2)),
       sigma = runif(1,5,20))
})

# Compile and fit model to stan
stanfit_Fabens <- stan(file = "Bayesian_Fabens.stan",
                data = data_listF,
                init = init_listF,
                chains = 4,
                iter = 50000,
                warmup = 25000,
                thin = 10)

# Results and diagnostics for Fabens growth model ----
stanfit_Fabens

# View visual diagnostics for the Fabens model
traceplot(stanfit_Fabens, pars = c("Linf", "K", "sigma")) +
  ggtitle("Traceplot of the Fabens growth model")

stan_dens(stanfit_Fabens, pars = c("Linf", "K", "sigma"))

# Get Credible Intervals for Linf and K
post_Fabens <- as.data.frame(stanfit_Fabens) %>% 
  summarise(Lci_lo = quantile(Linf, probs = .05),
            Lci_hi = quantile(Linf, probs = .95),
            Kci_lo = quantile(K, probs = .05),
            Kci_hi = quantile(K, probs = .95))
