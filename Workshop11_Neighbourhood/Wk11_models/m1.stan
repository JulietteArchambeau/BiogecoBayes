data{  
  int N ;                       // number of individuals  
  vector<lower=0> [N] dbh ;     // diameter at breast height
  vector [N] growth ;           // growth
  vector [N] D[N] ;             // distance matrix
  vector [N] S[N] ;             // size matrix
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
  vector[N] smat[N] ;
  vector[N] dmat[N] ;
  for(i in 1:N){
    for(j in 1:N){
     smat[i,j] = S[i,j]^gamma ;  
     dmat[i,j] = D[i,j]*D[i,j]*delta ;
    }
  }
  for(n in 1:N)
      NCI[n] = sum(smat[n] ./ exp(dmat[n])) ; // vectorial division
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
