data {                                                                         // observed variables 
  int<lower=1> N;                                                              // Number of observations
  vector[N] y;                                                                 // Response variable
  vector[N] age;                                                               // Tree age
  int<lower=0> nprov;                                                          // Number of provenances
  int<lower=0> nblock;                                                         // Number of blocks
  int<lower=0, upper=nprov> prov[N];                                           // Provenances
  int<lower=0, upper=nblock> bloc[N];                                          // Blocks
}



parameters {                                                            // unobserved variables
  real beta_age;
  real beta_age2;
  real alpha;
  
  real<lower=0> sigma_y;
  real<lower=0> sigma_block;                                            // stddev among intercepts of the blocks 
  vector[nblock] z_alpha_block;                                           // average intercept of each block
  
  corr_matrix[2] Rho_prov;
  vector<lower=0>[2] sigma_prov;
  vector[nprov] z_alpha_prov;
  vector[nprov] z_beta_prov;
}

transformed parameters {
  vector[2] v_prov[nprov];
  
  for(j in 1:nprov) v_prov[j] = [z_alpha_prov[j], z_beta_prov[j]]';
}

model{
  real mu[N];

//Priors
  Rho_prov ~ lkj_corr(4);
  sigma_prov ~ exponential(1);
  
  target += multi_normal_lpdf(v_prov|rep_vector(0, 2), Rho_prov);
  
  alpha ~ lognormal(0,1);
  beta_age ~ lognormal(0,1);
  beta_age2 ~ normal(0,1);
  
  z_alpha_block ~ normal(0,1);
  sigma_block ~ exponential(1);
  sigma_y ~ exponential(1);

// Linear model
  for (i in 1:N){
  mu[i] = alpha  + z_alpha_block[bloc[i]]*sigma_block + z_alpha_prov[prov[i]]*sigma_prov[1] + z_beta_prov[prov[i]]* sigma_prov[2] * age[i] + beta_age * age[i] + beta_age2 * age[i] * age[i];
  }  
  
// Likelihood
  y ~ lognormal(mu, sigma_y);                                           // or target += lognormal_lpdf(y|mu,sigma_y)
}

generated quantities {
  vector[N] y_rep;

  for(i in 1:N)  y_rep[i] = lognormal_rng(alpha  + z_alpha_block[bloc[i]]*sigma_block + z_alpha_prov[prov[i]]*sigma_prov[1] + z_beta_prov[prov[i]]* sigma_prov[2] * age[i] + beta_age * age[i] + beta_age2 * age[i] * age[i], sigma_y);
}

