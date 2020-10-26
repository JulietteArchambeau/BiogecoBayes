
data {                                             // observations
  int<lower=1> T;                                  // number of observations
  int<lower=0,upper=1> x[T];                       // observed variable
 
}
parameters {                                       // unobserved parameters
  simplex[2] Gamma[2];
  simplex[2] E[2]; 
}
model {
  //useful variables
  real acc[2];
  vector[2] alpha[T]; // T is the first dimension of the array here
  alpha[1,1:2] = rep_vector(0.5,2);

  //priors -- note that rows of the matrices must sum to one. Dirichlet priors. 
   for (n in 1:2){
   Gamma[n] ~ dirichlet(rep_vector(10, 2));
   E[n] ~ dirichlet(rep_vector(10, 2));
   }
  // hoping these will have rows that sum to 1 -- they should. 


  for (t in 2:T){ //loop over time
  
  for (k in 1:2){ //loop over categories for present state k
     for (j in 1:2){//loop over categories for past state j
           acc[j] = alpha[t - 1, j] * Gamma[j,k] * E[k,1+x[t]]; 
      }
      alpha[t,k] = sum(acc);
    }
  }

  target += log(sum(alpha[T])); //should there have been something more? Like, initial probabilities to estimate too? 
}

