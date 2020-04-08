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
  real<lower=0> beta_age;
  real beta_age2;
  vector<lower=0>[nprov] alpha_prov;                                           // Do not forget to define alpha_prov and alpha_block on R+ !
  vector<lower=0>[nblock] alpha_block;
  real<lower=0> sigma_y;
}


model{
  real mu[N];
  
  //Priors
  beta_age ~ lognormal(log(0.5), 1);
  beta_age2 ~ normal(0, 1);
  alpha_prov ~ normal(0, 1);
  alpha_block ~ normal(0, 1);
  sigma_y ~ exponential(1);
  
  // Likelihood
  for (i in 1:N){
    mu[i] = alpha_prov[prov[i]] + alpha_block[bloc[i]] + beta_age * age[i] + beta_age2 * square(age[i]);
  }
  y ~ lognormal(log(mu), sigma_y);                                            // mu is the median of the lognormal distribution
}

generated quantities {
  vector[N] y_rep;
  
  for(i in 1:N)  y_rep[i] = lognormal_rng(log(alpha_prov[prov[i]] + alpha_block[bloc[i]] + beta_age * age[i] + beta_age2 * square(age[i])), sigma_y);
}

