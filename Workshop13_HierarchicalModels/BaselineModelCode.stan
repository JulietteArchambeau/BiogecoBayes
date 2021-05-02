data {
  int<lower = 1> n_obs;                                              // Total number of observations
  int<lower = 1> n_age;                                              // Number of different age classes
  int<lower = 1> n_site;                                             // Number of different sites
  vector[n_obs] FLOWERING;                                           // Response variable (flowering dates)
  int<lower = 1, upper = n_age> AGE[n_obs];                          // Age variable
  int<lower = 1, upper = n_site> SITE[n_obs];                        // Site variable
  real prior_location_beta0;
  real<lower = 0.0> prior_scale_beta0;
}

parameters {
  simplex[3] pi;                                                     // unit complex specifying that the sum of its elements equal to one. 
  real beta0;                                                        // global intercept
  real<lower = 0.0> sigma_tot;                                       // Total standard deviation
  vector[n_age] alpha;                                               // Age intercepts
  vector[n_site] delta;                                              // Site intercepts
}

transformed parameters {
  real sigma;                                                        // Residual standard deviation
  real sigma_age;                                                    // Standard deviation of the age intercepts
  real sigma_site;                                                   // Standard deviation of the site intercepts
  vector[n_obs] mu;                                                  // linear predictor
  
  sigma = sqrt(pi[1]) * sigma_tot;
  sigma_age = sqrt(pi[2]) * sigma_tot;
  sigma_site = sqrt(pi[3]) * sigma_tot;
  mu = rep_vector(beta0, n_obs) + alpha[AGE] + delta[SITE];
    //for (i in 1:n_obs) {
    // mu[i] = beta0 + alpha[AGE[i]] + delta[SITE[i]];
    //}
}

model {
  // Priors
  beta0 ~ normal(prior_location_beta0, prior_scale_beta0);           // Prior of the global intercept
  sigma_tot ~ student_t(3, 0.0, 1.0);                                // Prior of the total standard deviation
  alpha ~ normal(0.0, sigma_age);                                    // Prior of the age intercepts
  delta ~ normal(0.0, sigma_site);                                   // Prior of the site intercepts
  
  // Likelihood
  FLOWERING ~ normal(mu, sigma);
}

generated quantities {
  vector[n_obs] log_lik;                                             // Log-likelihood
  vector[n_obs] y_rep;                                               // posterior predictive check

  for(i in 1:n_obs) {
    log_lik[i] = normal_lpdf(FLOWERING[i]| mu[i], sigma);
    y_rep[i] = normal_rng(mu[i], sigma);
  }
}
