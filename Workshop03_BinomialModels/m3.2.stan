data {                                             // observed variables 
  int<lower=1> N;                                  // number of observations
  int<lower=0,upper=50> y[N];                      // response variable
  int<lower=0,upper=50> z[N];                      // sample size per group  
  int<lower=0,upper=1>  A[N];
  int<lower=0,upper=1>  P[N];
  int<lower=0,upper=1>  V[N];
  
}

parameters {                                       // unobserved variables
  real alpha;
  real beta_P; 
  real beta_V; 
  real beta_A; 
  real gamma_b; 
}

model {
  for (k in 1:N){
  y[k] ~ binomial(z[k], inv_logit(alpha + beta_P*P[k] + beta_V*V[k] +beta_A*A[k]) );       // likelihood
  }
  alpha ~ normal(0, 1);                // prior of the mean 
  beta_P ~ normal(0, gamma_b);          // prior of the "slope" (here the effect, rather)
  beta_V ~ normal(0, gamma_b);          // prior of the slope
  beta_A ~ normal(0, gamma_b);          // prior of the slope
  gamma_b ~ gamma(0.1,0.1); 		// hyperprior
}
