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
  
  vector[nblock] z_alpha_block;                                                // z-score of block intercepts
  real<lower=0> sigma_block;                                                   // sd among intercepts of the blocks 
  
  vector<lower=0>[2] sigma_prov;                                               // sd among intercepts and slopes of the provenances
  matrix[2,nprov] z_prov;
  cholesky_factor_corr[2] L_R_prov;                                            // Cholesky correlation factors
}

transformed parameters {
  matrix[nprov,2] v_prov;
  vector[nprov] alpha_prov;
  vector[nprov] beta_prov;
  matrix[2,2] R_prov;
  
  v_prov = (diag_pre_multiply(sigma_prov,L_R_prov) * z_prov)';  
  alpha_prov  = v_prov[,1];
  beta_prov = v_prov[,2];
  R_prov = L_R_prov * L_R_prov';
}

model{
  real mu[N];

//Priors
  L_R_prov ~ lkj_corr_cholesky_lpdf(4);
  to_vector(z_prov) ~ normal_lpdf(0,1);
  sigma_prov ~ exponential(1);

  alpha ~ normal(0,1);
  beta_age ~ normal(0,1);
  beta_age2 ~ normal(0,1);
  
  z_alpha_block ~ normal(0,1);
  sigma_block ~ exponential(1);
  sigma_y ~ exponential(1);

// Linear model
  for (i in 1:N){
  mu[i] = alpha  + z_alpha_block[bloc[i]]*sigma_block + alpha_prov[prov[i]] + beta_prov[prov[i]] * age[i] + beta_age*age[i] + beta_age2*square(age)[i];
  }  
  
// Likelihood
  y ~ normal(mu, sigma_y);
}

generated quantities {
 
  vector[N] y_rep;
  for(i in 1:N)  y_rep[i] = normal_rng(alpha  + z_alpha_block[bloc[i]]*sigma_block + alpha_prov[prov[i]] + beta_prov[prov[i]] * age[i] + beta_age*age[i] + beta_age2*square(age)[i], sigma_y);
}

