data {                                             // observed variables 
  int<lower=1> N;                                  // number of observations
  int<lower=0,upper=50> y[N];                      // response variable
  int<lower=0,upper=50> z[N];                      // sample size per group  
  int<lower=0,upper=1>  A[N];                      // attribute age of pirating eagle
  int<lower=0,upper=1>  P[N];                      // attribute size of pirating eagle
  int<lower=0,upper=1>  V[N];                      // attribute size of victim eagle
  
}

parameters {                                       // unobserved variables
  real alpha;
  real beta_P; 
  real beta_V; 
  real beta_A; 
}

model {
  for (k in 1:N){
  y[k] ~ binomial(z[k], inv_logit(alpha + beta_P*P[k] + beta_V*V[k] +beta_A*A[k]) );       // likelihood
  }
  alpha ~ normal(0, 10);        // prior of the mean temp
  beta_P ~ normal(0, 5);          // prior of the slope
  beta_V ~ normal(0, 5);          // prior of the slope
  beta_A ~ normal(0, 5);          // prior of the slope
}

