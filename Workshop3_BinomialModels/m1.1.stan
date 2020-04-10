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
  y ~ binomial_logit(z, gamma*(temp - mu_temp)) ;
  mu_temp ~ normal(2, 10);        // prior of the mean temp
  gamma ~ normal(1, 10);          // prior of the slope
}
