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
  
  vector[nprov] z_alpha_prov;                                                  // z-score of prov intercepts
  vector[nprov] z_beta_prov;                                                   // z-score of prov slopes (with age)
  corr_matrix[2] R_prov;                                                       // correlation matrix R                               
  vector<lower=0>[2] sigma_prov;                                               // sd among intercepts and slopes of the provenances
}

transformed parameters {
  vector[2] v_prov[nprov];
  vector[nblock] alpha_block;
  vector[nprov] alpha_prov;
  vector[nprov] beta_prov;
  
  alpha_block =  z_alpha_block*sigma_block;
  alpha_prov = z_alpha_prov*sigma_prov[1];
  beta_prov = z_beta_prov*sigma_prov[2];
  
  for(j in 1:nprov) v_prov[j] = [z_alpha_prov[j], z_beta_prov[j]]';
}

model{
  real mu[N];

//Priors
  R_prov ~ lkj_corr(4);
  sigma_prov ~ exponential(1);
  
  target += multi_normal_lpdf(v_prov|rep_vector(0, 2), R_prov);
  
  alpha ~ normal(0,1);
  beta_age ~ normal(0,1);
  beta_age2 ~ normal(0,1);
  
  z_alpha_block ~ normal(0,1);
  sigma_block ~ exponential(1);
  sigma_y ~ exponential(1);

// Linear model
  for (i in 1:N){
  mu[i] = alpha  + alpha_block[bloc[i]] + alpha_prov[prov[i]] + beta_prov[prov[i]]* age[i] + beta_age*age[i] + beta_age2*square(age)[i];
  }  
  
// Likelihood
  y ~ normal(mu, sigma_y);
}

generated quantities {
  vector[N] y_rep;

  for(i in 1:N)  y_rep[i] = normal_rng(alpha  + alpha_block[bloc[i]] + alpha_prov[prov[i]] + beta_prov[prov[i]]* age[i] + beta_age*age[i] + beta_age2*square(age)[i], sigma_y);
}

