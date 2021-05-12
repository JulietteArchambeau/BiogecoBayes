functions {
  real TwoGaussianMixture_lpdf(real y, real prob, vector location, real scale) {
    real log_pdf[2];
    log_pdf[1] = log1m(prob) + normal_lpdf(y| location[1], scale);
    log_pdf[2] = log(prob) + normal_lpdf(y| location[2], scale);
    return log_sum_exp(log_pdf);
  }
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
  int<lower = 1> n_cultivar;					     // Number of different cultivars
  vector[n_obs] FLOWERING;                                           // Response variable (flowering dates)
  int<lower = 1, upper = n_age> AGE[n_obs];                          // Age variable
  int<lower = 1, upper = n_site> SITE[n_obs];                        // Site variable
  int<lower = 1, upper = n_cultivar> CULTIVAR[n_obs];		     // Cultivar variable
  real prior_location_beta0;
  real<lower = 0.0> prior_scale_beta0;
  real prior_location_diff;
  real<lower = 0.0> prior_scale_diff;
  real prior_location_eta0;
  real<lower = 0.0> prior_scale_eta0;
}

parameters {
  simplex[3] pi;                                                     // unit complex specifying that the sum of its elements equal to one. 
  real beta0;                                                        // global intercept
  real<lower = 0.0> sigma_tot;                                       // Total standard deviation
  vector[n_age] alpha;                                               // Age intercepts
  vector[n_site] delta;                                              // Site intercepts
  real<lower = 0> difference;                                        // difference between beta_1 and beta_2
  real<lower = 0.0> sigma_cultivar;
  real eta0;
  vector[n_cultivar] eta;
}

transformed parameters {
  real sigma;                                                        // Residual standard deviation
  real sigma_age;                                                    // Standard deviation of the age intercepts
  real sigma_site;                                                   // Standard deviation of the site intercepts
  vector[n_obs] mu[2];                                               // linear predictor
  vector[n_obs] p;						     // proba (early or late flowering)
  vector[2] beta;
  beta[1] = prior_location_beta0 + beta0 * prior_scale_beta0;
  beta[2] = beta[1] + difference;
  sigma = sqrt(pi[1]) * sigma_tot;
  sigma_age = sqrt(pi[2]) * sigma_tot;
  sigma_site = sqrt(pi[3]) * sigma_tot;
  mu[1] = beta[1] + alpha[AGE] + delta[SITE];
  mu[2] = beta[2] + alpha[AGE] + delta[SITE];
  p = inv_logit(rep_vector(eta0, n_obs) + eta[CULTIVAR]);
}

model {
  // Priors
  sigma_tot ~ student_t(3, 0.0, 1.0);
  sigma_cultivar ~ student_t(3, 0.0, 1.0);
  beta0 ~ normal(0.0, 1.0);
  difference ~ normal(prior_location_diff, prior_scale_diff);
  alpha ~ normal(0.0, sigma_age);
  delta ~ normal(0.0, sigma_site);
  eta0 ~ normal(prior_location_eta0, prior_scale_eta0);
  eta ~ normal(0.0, sigma_cultivar);

  // Likelihood
  for(i in 1:n_obs) {
    target += log_sum_exp(log1m(p[i]) + normal_lpdf(FLOWERING[i] | mu[1, i], sigma),
    log(p[i]) + normal_lpdf(FLOWERING[i] | mu[2, i], sigma));
  }
}

generated quantities {
  vector[n_obs] log_lik;                                             // Log-likelihood
  vector[n_obs] y_rep;                                               // posterior predictive check

  for(i in 1:n_obs) {
    log_lik[i] = log_sum_exp(log1m(p[i]) + normal_lpdf(FLOWERING[i] | mu[1, i], sigma),
    log(p[i]) + normal_lpdf(FLOWERING[i] | mu[2, i], sigma));
    y_rep[i] = TwoGaussianMixture_rng(p[i], to_vector(mu[1:2, i]), sigma);
  }
}
