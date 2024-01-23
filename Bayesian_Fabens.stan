// Bayesian version of the growth model after Fabens (1965)
data {
  int<lower=0> N;  // number of observations
  real Lm[N];      // length at mark
  real Lr[N];      // length at recapture
  real dt[N];      // time difference between mark and recapture
  real mLinf;      // prior for mean Linf on log scale
  real mK;         // prior for mean K on log scale
  real sdLinf;     // prior for standard deviation of Linf on a log scale
  real sdK;        // prior for standard deviation of K on a log scale
}

parameters {
  real<lower=0> sigma;                  // additive Gaussian error
  real<lower=-10, upper=10> logLinf;    // log scale Linf
  real<lower=-5, upper=5> logK;         // log scale K
}

transformed parameters {
  real Linf;
  real K;
  
  Linf = exp(logLinf);                  // transforming parameter back to original scale 
  K = exp(logK);                        // transforming parameter back to original scale 
}

model {
  vector[N] ypred;
  
  sigma ~ uniform(0, 100);              // reference prior for the standard deviation
  logLinf ~ normal(mLinf, sdLinf);      // prior probability distribution for Linf
  logK ~ normal(mK, sdK);               // prior probability distribution for K
  
  // calculate likelihood for data
  for(i in 1:N){
    ypred[i] = Lm[i] + (Linf-Lm[i]) * (1-exp(-K*dt[i]));
  }
  Lr~normal(ypred, sigma);
}

generated quantities { // predicted values for inspecting model fit
  vector[N] predy;
  for(i in 1:N){
    predy[i] = normal_rng(Lm[i] + (Linf-Lm[i]) * (1-exp(-K*dt[i])), sigma);
  }
}
