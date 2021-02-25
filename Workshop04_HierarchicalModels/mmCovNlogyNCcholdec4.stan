data {
  int<lower=1> N ;
  vector[N] y ;
  vector[N] age ;
  int<lower=1> nprov ;
  int<lower=1> nblock ;
  int<lower=1, upper=nprov> prov[N] ;
  int<lower=1, upper=nblock> bloc[N] ;
}
parameters {
  real alpha ;
  real beta[2] ;
  vector[nblock] z_alpha_block ;
  matrix[2,nprov] z_prov ;
  vector<lower=0>[4] sigma ;
  cholesky_factor_corr[2] L_R_prov ;
}
transformed parameters {
  matrix[nprov,2] v_prov = (diag_pre_multiply(sigma[3:4],L_R_prov) * z_prov)' ;
  vector[nprov] alpha_prov = v_prov[,1] ;
  vector[nprov] beta_prov = v_prov[,2] ;
  matrix[2,2] R_prov = L_R_prov * L_R_prov' ;
}
model{
  y ~ normal(alpha + 
               z_alpha_block[bloc]*sigma[2] + 
               alpha_prov[prov] + 
               beta_prov[prov] .* age + 
               beta[1]*age + 
               beta[2]*square(age),
             sigma[1]) ;
  to_vector(z_prov) ~ std_normal() ;
  z_alpha_block ~ std_normal() ;
  alpha ~ std_normal() ;
  beta ~ std_normal() ;
  sigma ~ exponential(1) ;
  R_prov ~ lkj_corr(4) ;
}
generated quantities {
  vector[N] y_rep = alpha + 
                      z_alpha_block[bloc]*sigma[2] + 
                      alpha_prov[prov] + 
                      beta_prov[prov] .* age + 
                      beta[1]*age + 
                      beta[2]*square(age) ;
}

