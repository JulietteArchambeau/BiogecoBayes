data{  
  int N ;                       // number of individuals  
  int M ;                       // number of distances
  vector<lower=0> [N] dbh ;     // diameter at breast height
  vector<lower=0> [N] growth ;  // growth
  vector<lower=0> [M] d ;       // distances
  int nb_b[N] ;                 // non zero values
  int pos[N] ;                  // first non zero position
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
                                      - delta*to_vector(segment(d, pos[n], nb_b[n])))),
                    sigma) ;
}
