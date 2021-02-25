functions {
  real dirichlet_multinomial_lpmf(int[] y, vector alpha) {
    real alpha_plus = sum(alpha);
    return lgamma(alpha_plus) + sum(lgamma(alpha + to_vector(y)))
                - lgamma(alpha_plus+sum(y)) - sum(lgamma(alpha));
  }
}
data {
  int<lower = 1> N  ;       // # of individuals
  int<lower = 1> S  ;       // # of species
  int<lower = 1> K  ;       // # of environmental descriptors
  int<lower = 0, upper=1> Y[N, S]  ; // individuals presence or absence for each species
  matrix[N,K] X ;           // environmental descriptors
}
parameters {
  vector[S] alpha ; // intercept
  matrix[S,K] beta ; // sigmoidal slope
  matrix[S,K] gamma ; // quadratic form
}
model {
  alpha ~ normal(0,1) ;  
  for(k in 1:K){
    beta[,k] ~ normal(0,1) ;
    gamma[,k] ~ normal(0,1) ;
  }
  for (n in 1:N)
    Y[n] ~ dirichlet_multinomial(softmax(alpha + beta*to_vector(X[n,]) + gamma*to_vector(X[n,] .* X[n,])))  ; // likelihood
}
