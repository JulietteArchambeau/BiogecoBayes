functions {
  real TwoGaussianMixture_rng(real prob, vector location, real scale) {
    int z;
    z = bernoulli_rng(prob);
    return  z ? normal_rng(location[2], scale) : normal_rng(location[1], scale); 
  }
}

data {
  int<lower = 1> n_obs;                                              // Total number of observations
  int<lower = 1> n_age;                                              // Number of different age classes
  int<lower = 1> n_site;                                             // Number of different sites
  vector[n_obs] FLOWERING;                                           // Response variable (flowering dates)
  int<lower = 1, upper = n_age> AGE[n_obs];                          // Age variable
  int<lower = 1, upper = n_site> SITE[n_obs];                        // Site variable
  real prior_location_beta0;
  real<lower = 0.0> prior_scale_beta0;
  real prior_location_diff;
  real<lower = 0.0> prior_scale_diff;
}

parameters {
  real<lower = 0.0, upper = 1.0> p;				     // proba (early or late flowering)
  simplex[3] pi;                                                     // unit complex specifying that the sum of its elements equal to one. 
  real beta0;                                                        // global intercept
  real<lower = 0.0> sigma_tot;                                       // Total standard deviation
  vector[n_age] alpha;                                               // Age intercepts
  vector[n_site] delta;                                              // Site intercepts
  real<lower = 0> difference;                                        // difference between beta_1 and beta_2
}

transformed parameters {
  real sigma;                                                        // Residual standard deviation
  real sigma_age;                                                    // Standard deviation of the age intercepts
  real sigma_site;                                                   // Standard deviation of the site intercepts
  vector[n_obs] mu[2];                                               // linear predictor
  vector[2] beta;
  beta[1] = prior_location_beta0 + beta0 * prior_scale_beta0;
  beta[2] = beta[1] + difference;
  sigma = sqrt(pi[1]) * sigma_tot;
  sigma_age = sqrt(pi[2]) * sigma_tot;
  sigma_site = sqrt(pi[3]) * sigma_tot;
  mu[1] = rep_vector(beta[1], n_obs) + alpha[AGE] + delta[SITE];
  mu[2] = rep_vector(beta[2], n_obs) + alpha[AGE] + delta[SITE];
  //for (i in 1:n_obs) {
  // mu[1, i] = beta[1] + alpha[AGE[i]] + delta[ZONE[i]];
  // mu[2, i] = beta[2] + alpha[AGE[i]] + delta[ZONE[i]];
  //}
}

model {
  // Priors
  beta0 ~ normal(0.0, 1.0);
  difference ~ normal(prior_location_diff, prior_scale_diff);
  sigma_tot ~ student_t(3, 0.0, 1.0);                                // Prior of the total standard deviation
  alpha ~ normal(0.0, sigma_age);                                    // Prior of the age intercepts
  delta ~ normal(0.0, sigma_site);                                   // Prior of the site intercepts
  
  // Likelihood
  for(i in 1:n_obs) {
    target += log_sum_exp(log1m(p) + normal_lpdf(FLOWERING[i] | mu[1, i], sigma), log(p) + normal_lpdf(FLOWERING[i] | mu[2, i], sigma));
  }
}

generated quantities {
  vector[n_obs] log_lik;                                             // Log-likelihood
  vector[n_obs] y_rep;                                               // posterior predictive check

  for(i in 1:n_obs) {
    log_lik[i] = log_sum_exp(log1m(p) + normal_lpdf(FLOWERING[i] | mu[1, i], sigma), log(p) + normal_lpdf(FLOWERING[i] | mu[2, i], sigma));
    y_rep[i] = TwoGaussianMixture_rng(p, to_vector(mu[1:2, i]), sigma);
  }
}
