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
  matrix[2,nprov] z_prov;                                                      // z-score of the provenance intercepts and slopes
  cholesky_factor_corr[2] L_R_prov;                                            // Cholesky correlation factors
}

transformed parameters {
  matrix[2,nprov] v_prov;
  v_prov = (diag_pre_multiply(sigma_prov,L_R_prov) * z_prov);   // Here, we do not transpose v_prov!
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
mu[i] = alpha  + z_alpha_block[bloc[i]]*sigma_block + v_prov[1,prov[i]] + v_prov[2,prov[i]] * age[i] + beta_age*age[i] + beta_age2*square(age)[i];
// Careful here: v_prov[effect,prov] instead of v_prov[prov,effect]
}  

// Likelihood
y ~ normal(mu, sigma_y);
}

generated quantities {

// Return the correlation matrix
matrix[2,2] R_prov;
vector[N] y_rep;

R_prov = multiply_lower_tri_self_transpose(L_R_prov); // transform the Cholesky factor into an ordinary correlation matrix

// Return simulated data
for(i in 1:N) { 
y_rep[i] = normal_rng(alpha  + z_alpha_block[bloc[i]]*sigma_block + v_prov[1,prov[i]] + v_prov[2,prov[i]] * age[i] + beta_age*age[i] + beta_age2*square(age)[i], sigma_y);
}
}

