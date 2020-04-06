data {                                                                         // observed variables 
  int<lower=1> N;                                                              // Number of observations
  vector[N] y;                                                                 // Response variable
  vector[N] age;                                                               // Tree age
  int<lower=0> nprov;                                                          // Number of provenances
  int<lower=0> nblock;                                                         // Number of blocks
  int<lower=0, upper=nprov> prov[N];                                           // Provenances
  int<lower=0, upper=nblock> bloc[N];                                          // Blocks
}



parameters {                                                                   // unobserved variables
  real beta_age;
  real beta_age2;
  real alpha;
  
  real<lower=0> sigma_y;
  real<lower=0> sigma_block;                                             // stddev among intercepts of the blocks 
  vector[nblock] alpha_block;                                                  // average intercept of each block
  
  matrix[2,nprov] z_u;
  vector<lower=0>[2] sigma_u;
  cholesky_factor_corr[2] L_u;
  
}

transformed parameters {
  matrix[2,nprov] u;
  
  u = diag_pre_multiply(sigma_u,L_u)*z_u;
}

model{
  real mu[N];

//Priors
  L_u ~ lkj_corr_cholesky(2);
  to_vector(z_u) ~ normal(0,1);
  
  alpha ~ normal(0,1);
  beta_age ~ normal(0,1);
  beta_age2 ~ normal(0,1);
  
  alpha_block ~ normal(0, sigma_block);
  sigma_block ~ cauchy(0,1);
  sigma_y ~ cauchy(0,1);

// Likelihood
  for (i in 1:N){
  mu[i] = alpha  + alpha_block[bloc[i]] + u[1,prov[i]] + u[2,prov[i]] * age[i] + beta_age * age[i] + beta_age2 * age[i] * age[i];
  }
  y ~ lognormal(mu, sigma_y);
}

generated quantities {
  vector[N] y_rep;

  for(i in 1:N)  y_rep[i] = lognormal_rng(alpha  + alpha_block[bloc[i]] + u[1,prov[i]] + u[2,prov[i]] * age[i] + beta_age * age[i] + beta_age2 * age[i] * age[i], sigma_y);
}


