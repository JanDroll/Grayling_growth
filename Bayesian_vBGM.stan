// Bayesian version of the von Bertalanffy growth model
data {
  int<lower=0> N;
  real TL[N];
  int Age[N];
  real mLinf;
  real sdLinf;
  real mK;
  real sdK;
  real mt0;
  real sdt0;
}

parameters {
  real<lower=0> sigma;
  real<lower=-10, upper=10> logLinf;
  real<lower=-5, upper=5> logK;
  real<lower=-5, upper=5> t0;
}

transformed parameters {
  real Linf;
  real K;
  
  Linf = exp(logLinf);
  K = exp(logK);
}

model {
  vector[N] ypred;
  
  sigma ~ uniform(0, 100);
  logLinf ~ normal(mLinf, sdLinf);
  logK ~ normal(mK, sdK);
  t0 ~ normal(mt0, sdt0);
  
  for(i in 1:N){
    ypred[i] = Linf * (1-exp(-K * (Age[i]-t0)));
  }
}

generated quantities {
  vector[N] predy;
  for(i in 1:N){
    predy[i] = normal_rng(Linf * (1-exp(-K * (Age[i]-t0))), sigma);
  }
}
