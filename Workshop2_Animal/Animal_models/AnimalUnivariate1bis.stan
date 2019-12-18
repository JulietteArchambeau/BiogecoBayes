data {
  int<lower=0>  N ; // # of individuals
  real Y[N] ; // phenotype
  cov_matrix[N] K ; // kinship covariance matrix
}
transformed data{
  matrix[N, N] A ; // cholesky-decomposed kinship
  real<lower=0> sigma ; // phenotypic variance
  A = cholesky_decompose(K) ;
  sigma = sd(Y) * sd(Y) ;
}
parameters {
  vector[N]  u_tilde ; // random effects / breeding values
  real mu ; // intercept
  simplex[2] part ; // variance partition between environement and genetic
}
model {
    u_tilde ~ normal(0, 1) ; // priors
    mu ~ normal(0, 1) ;
    Y ~ normal(mu + sqrt(sigma*part[1])*(A * u_tilde), sqrt(sigma*part[2]));
}
generated quantities{
  real sigmaE ; // environmental variation
  real sigmaG ; // genetic variation
  sigmaE = sigma*part[2] ;
  sigmaG = sigma*part[1] ;
}
