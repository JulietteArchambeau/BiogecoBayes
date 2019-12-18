data {
  int<lower=0>  N ; // # of individuals
  int<lower=0>  J ; // # of covariates + 1 (intercept)  
  real Y[N] ; // phenotype
  matrix[N,J] X; // covariates
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
  vector[J] beta; // fixed effects
  simplex[2] part ; // variance partition between environement and genetic
}
model {
    u_tilde ~ normal(0, 1) ; // priors
    to_vector(beta) ~ normal(0, 1) ;
    Y ~ normal(X*beta + sqrt(sigma*part[1])*(A * u_tilde), sqrt(sigma*part[2])) ;
}
generated quantities{
  real sigmaE ; // environmental variation
  real sigmaG ; // genetic variation
  sigmaE = sigma*part[2] ;
  sigmaG = sigma*part[1] ;
}
