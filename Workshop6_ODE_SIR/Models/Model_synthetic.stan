//This bloc as been taken from jrmihalj.github.io/estimating-transmission-by-fitting-mechanistic-models-in-Stan/
functions {
  // This bloc is specific to ODE estimate : this function must contain the time derivative, and is very similar to the function we defined for the deSolve package.
  // This largely follows the deSolve package, but also includes the x_r and x_i variables.
  // These variables are mandatory, due to the signature of the derivative functions in STAN ODE integrators. They may contain for exemple additional parameters that we do not want to infer, and are kept fixed. Here, they will remain null (see transformed_data bloc)
  
   real[] SI(real t,
            real[] y,
            real[] params,
            real[] x_r,
            int[] x_i) {
      
      real dydt[3];
      
      dydt[1] = - params[1] * y[1] * y[2];
      dydt[2] = params[1] * y[1] * y[2] - params[2] * y[2];
      dydt[3] = params[2] * y[2];
      
      return dydt;
    }
  
}

data {
  int<lower = 1> n_obs; // Number of days sampled
  int<lower = 1> n_params; // Number of model parameters
  int<lower = 1> n_difeq; // Number of differential equations in the system
  int<lower = 1> n_sample; // Number of hosts sampled at each time point.
  int<lower = 1> n_fake; // This is to generate "predicted"/"unsampled" data
  
  int y[n_obs]; // The binomially distributed data
  real t0; // Initial time point (zero)
  real ts[n_obs]; // Time points that were sampled
  
  real fake_ts[n_fake]; // Time points for "predicted"/"unsampled" data
}

transformed data {
  real x_r[0]; // No additional real parameter to ODE time derivative function => is set to a vector of dimension 0
  int x_i[0]; // No additional integer parameter to ODE time derivative function => is set to a vector of dimension 0
}

parameters {
  real<lower = 0> params[n_params]; // Model parameters
  real<lower = 0, upper = 1> S0; // Initial fraction of hosts susceptible
}

transformed parameters{
  real y_hat[n_obs, n_difeq]; // Output from the ODE solver
  real y0[n_difeq]; // Initial conditions for both S and I

  y0[1] = S0;
  y0[2] = 1 - S0;
  y0[3] = 0;
  
  y_hat = integrate_ode_rk45(SI, y0, t0, ts, params, x_r, x_i);
  
}

model {
  params ~ normal(0, 2); //constrained to be positive
  S0 ~ normal(0.5, 0.5); //constrained to be 0-1.
  
  y ~ binomial(n_sample, y_hat[, 2]); //y_hat[,2] are the fractions infected from the ODE solver
  
}

generated quantities {
  // Generate predicted data over the whole time series:
  real fake_I[n_fake, n_difeq];
  
  fake_I = integrate_ode_rk45(SI, y0, t0, fake_ts, params, x_r, x_i);
  
}




