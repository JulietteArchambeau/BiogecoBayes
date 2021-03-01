data{  
  int N ;                       // number of individuals  
  int M ;                       // number of distances
  vector<lower=0> [N] dbh ;     // diameter at breast height
  vector [N] growth ;           // growth
  vector [M] d ;                // distance matrix
  vector [M] s ;                // size matrix
  int n_nb[N] ;                 // non zero values
  int pos[N] ;                  // first non zero position
}
parameters{
  real alpha ;            // Intercept
  real beta[2] ;          // DBH and NCI slopes
  real<lower=0> gamma ;   // neighbour diameter effect
  real<lower=0> delta ;   // neighbour distance effect
  real<lower=0> sigma ;   // residual variation
}
transformed parameters{  
  vector[N] NCI ;  
  vector[M] svec ;
  vector[M] dvec ;
  for(m in 1:M){
     svec[m] = s[m]^gamma ;  
     dvec[m] = d[m]*d[m]*delta ;
  }
  for(n in 1:N)
    NCI[n]=sum(exp(log(segment(svec, pos[n], n_nb[n])) -
                  (segment(dvec, pos[n], n_nb[n])))) ;
}
model{
  alpha ~ normal(0, 1) ;
  beta ~ normal(0, 1) ;
  gamma ~ normal(0, 1) ; 
  delta ~ normal(0, 1) ; 
  sigma ~ exponential(1) ;
  growth ~ normal(alpha + 
                  beta[1] * dbh +
                  beta[2] * NCI,
                  sigma) ;
}

