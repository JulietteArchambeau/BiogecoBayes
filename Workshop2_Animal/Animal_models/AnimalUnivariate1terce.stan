data {
  int<lower=0>  N ; // # of individuals
  real Y[N] ; // phenotype
  cov_matrix[N] K ; // kinship covariance matrix
}
transformed data{
  real<lower=0> sigma ; // phenotypic variance
  sigma = sd(Y) * sd(Y) ;
}
parameters {
  vector[N]  u ; // random effects / breeding values
  real mu ; // intercept
  simplex[2] part ; // variance partition between environement and genetic
}
model {
    u ~ multi_normal(rep_vector(0, N), sqrt(sigma*part[1])*K) ;
    mu ~ normal(0, 1) ;
    Y ~ normal(mu + u, sqrt(sigma*part[2]));
}
