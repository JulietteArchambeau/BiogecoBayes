data{  
  int N ;                       // number of individuals  
  vector<lower=0> [N] dbh ;     // diameter at breast height
  vector [N] growth ;           // growth
  matrix<lower=0> [N, N] D ;    // distance matrix
}
parameters{
  real<lower=0> alpha ;   // Intercept  
  real beta[2] ;          // DBH and NCI slopes
  real<lower=0> delta ;   // neighbour distance effect
  real<lower=0> gamma ;   // neighbour diameter effect
  real<lower=0> sigma ;   // residual variation
}
model{
  for(n in 1:N)
    growth ~ normal(alpha + 
                    beta[1] * log(dbh) + 
                    beta[2] * sum(exp(gamma*log(dbh) 
                                      - delta*to_vector(D[n,]))),
                    sigma) ;
}
