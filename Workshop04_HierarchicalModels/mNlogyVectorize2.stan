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
  vector[nprov] alpha_prov;
  vector[nblock] alpha_block;
  real<lower=0> sigma_y;
}


model{
  
  vector[N] mu;
  
  //Priors
  target += normal_lpdf(alpha_block | 0, 10);
  target += normal_lpdf(alpha_prov | 0, 10);
  target += normal_lpdf(beta_age | 0, 10);
  target += normal_lpdf(beta_age2 | 0, 10);
  target += cauchy_lpdf(sigma_y | 0, 25);
  
  // Likelihood
  target +=  normal_lpdf(y |alpha_prov[prov] + alpha_block[bloc] + beta_age * age + beta_age *square(age), sigma_y);
}

generated quantities {
  real y_rep[N] = normal_rng(alpha_prov[prov] + alpha_block[bloc] + beta_age * age + beta_age2*square(age), sigma_y);
}
