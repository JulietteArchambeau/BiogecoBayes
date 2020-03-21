data {                                             // observed variables 
  int<lower=1> N;                                  // number of observations
  vector[N] temp;                                  // temperature
  int<lower=0,upper=50> y[N];                      // response variable
  int<lower=0,upper=50> z[N];                      // sample size per group  
}
parameters {                                       // unobserved variables
  real mu_temp;
  real gamma; 
}
model {
  for (k in 1:N){
  y[k] ~ binomial(z[k], inv_logit(gamma*(temp[k]-mu_temp)));       // likelihood
  //y[k] ~ bernoulli_logit(alpha+beta*temp[k]);       // likelihood -- does not work because not bernoulli of course!
  }
  mu_temp ~ normal(2, 10);        // prior of the mean temp
  gamma ~ normal(1, 10);          // prior of the slope
}