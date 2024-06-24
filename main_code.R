library(tidyverse)
library(FSA)
library(TropFishR)
library(fishboot)
library(car)
library(rstan)
library(tidybayes)
library(bayesplot)
library(flextable)

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
                       date_season = col_date(format = "%d.%m.%y")), trim_ws = TRUE) %>% 
  mutate(length_cm = length/10)

# Generating mean and sd for the prior probability distribution of Linf, K, and t0
# -> shared for typical vB and Fabens
mean_Linf <- mean(log(priors$Linf))
sd_Linf <- sd(log(priors$Linf))
mean_K <- mean(log(priors$K))
sd_K <- sd(log(priors$K))
mean_t0 <- mean(priors$t0)
sd_t0 <- sd(priors$t0)

## Frequentist von Bertalanffy Growth model ----
svTypical <- vbStarts(tl~age, data = AL_gray, plot = TRUE)
vb <- vbFuns("Typical")
f.fit <- nls(tl ~ vb(age, Linf, K, t0), data = AL_gray, start = svTypical)

# Getting CIs for the parameters of the von Bertalanffy growth model
coef(f.fit)
vbBoot <- Boot(f.fit)
confint(vbBoot, level = .95)

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
init_listVB <- lapply(1:11, function(i) {
  list(logLinf = log(rnorm(1, max(data_listVB$TL), 1)),
       logK = log(runif(1,0.01,5)),
       t0 = rnorm(1,1,1),
       sigma = runif(1,1,10))
})

# Compile nd fit model to stan
# Important! : If using this code on another computer the number of chains 
# must correspond to the number of cores of the CPU. The number of cores can 
# be found out by using parallel::detectCores(). The number of initial values
# has to be adjusted as well in the lapply function above.
stanfit_vBert <- stan(file = "Bayesian_vBGM.stan",
                      data = data_listVB,
                      init = init_listVB,
                      chains = 11,
                      iter = 50000,
                      warmup = 25000,
                      thin = 10)

# Results and diagnostics for Bayesian von Bertalanffy growth model ----
stanfit_vBert

# View visual diagnostics for the Fabens model
traceplot(stanfit_vBert, pars = c("Linf", "K", "t0", "sigma"), alpha = .7) +
  ggtitle("Traceplot of the von Bertalanffy growth model")

stan_dens(stanfit_vBert, pars = c("Linf", "K", "t0", "sigma")) +
  ggtitle("Density of the parameters for the Bayesian von Bertalanffy growth model")

pairs(stanfit_vBert, pars = c("Linf", "K", "t0", "sigma"))

# Get Credible Intervals for Linf, K, and t0
post_vB <- as.data.frame(stanfit_vBert) %>% 
  summarise(Lci_lo = quantile(Linf, probs = .05),
            Lci_hi = quantile(Linf, probs = .95),
            Kci_lo = quantile(K, probs = .05),
            Kci_hi = quantile(K, probs = .95),
            t0ci_lo = quantile(t0, probs = .05),
            t0ci_hi = quantile(t0, probs = .95))

## Frequentist Fabens growth model ----
# Finding starting values for Linf and K
max(Gray_mrt$Lr) # Starting value for Linf
with(Gray_mrt,mean((log(Lr)-log(Lm))/dt)) # Starting value for K

# Fitting the model with nls
Fabens.sv <- list(Linf = 504, K = 0.3945063)
Fabens_vb <- vbFuns("Fabens2")
FVB <- nls(Lr ~ Fabens_vb(Lm,dt,Linf,K), start = Fabens.sv, data = Gray_mrt)

# Bootstrapping for CIs
coef(FVB)
fvbBoot <- Boot(FVB)
confint(fvbBoot, level = .95)

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
init_listF <- lapply(1:11, function(i) {
  list(logLinf = log(rnorm(1, max(data_listF$Lr), 1)),
       logK = log(runif(1, mean(log(data_listF$Lr) - log(data_listF$Lm)) / data_listF$dt, 2)),
       sigma = runif(1,5,20))
})

# Compile and fit model to stan
# Important! : If using this code on another computer the number of chains 
# must correspond to the number of cores of the CPU. The number of cores can 
# be found out by using parallel::detectCores(). The number of initial values
# has to be adjusted as well in the lapply function above.
stanfit_Fabens <- stan(file = "Bayesian_Fabens.stan",
                data = data_listF,
                init = init_listF,
                chains = 11,
                iter = 50000,
                warmup = 25000,
                thin = 10)

# Results and diagnostics for Fabens growth model ----
stanfit_Fabens

# View visual diagnostics for the Fabens model
traceplot(stanfit_Fabens, pars = c("Linf", "K", "sigma"), ncol = 2) +
  ggtitle("Traceplot of the Fabens growth model")

stan_dens(stanfit_Fabens, pars = c("Linf", "K", "sigma"), ncol = 2) +
  ggtitle("Density of the parameters for the Bayesian Fabens growth model")

pairs(stanfit_Fabens, pars = c("Linf", "K", "sigma"))

# Get Credible Intervals for Linf and K
post_Fabens <- as.data.frame(stanfit_Fabens) %>% 
  summarise(Lci_lo = quantile(Linf, probs = .05),
            Lci_hi = quantile(Linf, probs = .95),
            Kci_lo = quantile(K, probs = .05),
            Kci_hi = quantile(K, probs = .95))

## ELEFAN with grayling data ----
# Creating lfq data object
lfq_new <- lfqCreate(data = lfq_data, Lname = "length_cm", Dname = "date_season")

plot(lfq_new, Fname = "catch", hist.sc = 1)

gray_MA <- lfqRestructure(lfq_new, MA = 7)
plot(gray_MA, hist.sc = 1)

# ELEFAN with genetic algorithm
gray_elefan.ga <- ELEFAN_GA(lfq_new,
                            seasonalised = FALSE,
                            low_par = list(Linf = 40, K = 0.1, t_anchor = 0, ts = 0, C = 0),
                            up_par = list(Linf = 60, K = 0.9, t_anchor = 1, ts = 1, C = 1),
                            popSize = 100,
                            maxiter = 100,
                            pmutation = 0.2,
                            run = 30,
                            MA = 7,
                            plot.score = TRUE,
                            plot = TRUE,
                            monitor = FALSE,
                            parallel = TRUE)

unlist(gray_elefan.ga$par)

## Results table ----
Growth_param <- tibble(Method = c("von Bertalanffy", "von Bertalanffy", "Fabens", "Fabens", "ELEFAN"),
                       Linf = c(519.56, 3785.10, 536.82, 538.85, 430.48),
                       L_up = c(586.85, 10137.64, 586.37, 678.30, 500.00),
                       L_low = c(460.88, 2021.51, 493.91, 474.32, 360.00),
                       K = c(0.21, 0.01, 0.38, 0.37, 0.12),
                       K_up = c(0.27, 0.03, 0.45, 0.49, 0.16),
                       K_low = c(0.17, 0.005, 0.31, 0.25, 0.10),
                       Type = c("Bayesian", "Frequentist", "Bayesian", "Frequentist", "Frequentist"),
                       Uncertainty = c("Credible Interval", "Confidence Interval", "Credible Interval", "Confidence Interval", "Confidence Interval"))

flextable(Growth_param) %>% 
  save_as_docx(path = "Growth_para.docx")
