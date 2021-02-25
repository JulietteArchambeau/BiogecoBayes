data {
  int<lower=1> N ;                // # of individuals
  int<lower=1> K ;                // # of environmental descriptors
  int<lower=0, upper=1> Y[N] ;    // individuals presence or absence (0-1) 
  matrix[N,K] X ;                 // environmental descriptors
  vector<lower=0, upper=1>[N] w ; // weights for pseudo-absence to be equal to presence
}
parameters {
  real alpha ;       // intercept
  vector[K] beta ;   // sigmoidal slope
  vector[K] gamma ;  // quadratic form
}
model {
  alpha ~ normal(0, 10^6) ; // priors
  beta ~ normal(0, 10^6) ;
  gamma ~ normal(0, 10^6) ;
  target += w*bernoulli_lpmf(Y | inv_logit(alpha + X * beta + X .* X * gamma)) ; // likelihood
}