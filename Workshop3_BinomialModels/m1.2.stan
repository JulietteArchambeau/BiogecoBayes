data {                                             // observed variables 
  int<lower=1> N;                                  // number of observations
  int<lower=0,upper=1> y[N];                      // response variable
}
parameters {                                       // unobserved variables
  real p;
  }
model {
  y ~ bernoulli(p);       // likelihood -- does not work because not bernoulli of course!
  p ~ beta(2,2);              // prior on the probability 
}