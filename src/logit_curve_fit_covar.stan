data {
  int<lower=1> N;          // Number of observations
  vector[N] x;             // Predictor variable (e.g., titre)
  vector[N] y;             // Continuous response variable
  array[N] int x_cov; // Covariate matrix (for hierarchical effect)
  int N_cov; // Number of levels

}
transformed data {
       
    real x_min = min(x);
    real x_max = max(x);
    real step = (x_max - x_min) / 99;
    real midpoint = (x_min + x_max) / 2;
}

parameters {
  real L_global;        // Global upper asymptote (shared prior)
  real<lower = 0> k_global;                  // Global steepness (shared prior)
  real x0_global;                          // Global midpoint (shared prior)
  real<lower=0> sigma;                     // Standard deviation of errors
  
  // Group-level deviations from the global mean
  array[N_cov] real L_deviation;  // Deviations for L per group
  array[N_cov] real k_deviation;  // Deviations for k per group
  array[N_cov] real x0_deviation; // Deviations for x0 per group
  
  // Standard deviations for the group-level deviations (hyperpriors)
  real<lower=0> L_sigma;  
  real<lower=0> k_sigma;  
  real<lower=0> x0_sigma;  
}

model {
  vector[N] y_hat;
  
  // Logistic function for curve fitting
  for (n in 1:N) {
    real L = inv_logit(L_global + L_deviation[x_cov[n]] * L_sigma);  // L for each group
    real k = k_global + k_deviation[x_cov[n]] * k_sigma;  // k for each group
    real x0 = x0_global + x0_deviation[x_cov[n]] * x0_sigma; // x0 for each group

    y_hat[n] = L * (1 - 1 / (1 + exp(-k * (x[n] - x0))));
  }

  // Likelihood: Assume normal residuals
  y ~ normal(y_hat, sigma);

  // Priors
  L_global ~ normal(0, 1.84);      // Prior for upper asymptote
  k_global ~ uniform(0, 5);        // Prior for steepness
  x0_global ~ normal(midpoint, midpoint / 2);       // Prior for midpoint

// Group-level priors (for the deviations)
  L_sigma ~ exponential(1);
  k_sigma ~ normal(0, 1);
  x0_sigma ~  exponential(1);

  // Deviations for each group drawn from normal distributions centered around 0
  L_deviation ~ normal(0, 1); // Priors for deviations of L
  k_deviation ~ normal(0, 1); // Priors for deviations of k
  x0_deviation ~ normal(0, 1); // Priors for deviations of x0
}
generated quantities {
    vector[100] y_hat_rel;
    vector[100] y_hat_new;
    vector[100] y_prot;

    matrix[100, N_cov] y_hat_rel_age;
    matrix[100, N_cov] y_hat_new_age;
    matrix[100, N_cov] y_prot_age;

    vector[100] x_new;

    vector[N_cov] L_temp;
    vector[N_cov] k_temp;
    vector[N_cov] x0_temp;
    
    // Posterior prediction
   
    for (j in 1:N_cov) {
      L_temp[j] = inv_logit(L_global + L_deviation[j] * L_sigma);  // L for each group
      k_temp[j] = k_global + k_deviation[j] * k_sigma;  // k for each group
      x0_temp[j] = x0_global + x0_deviation[j] * x0_sigma; // x0 for each group

      for (i in 1:100) {
          x_new[i] = x_min + (i - 1) * step;
          y_prot[i] =  1 / (1 + exp(-k_global * (x_new[i] - x0_global)));
          y_hat_new[i] = inv_logit(L_global) * (1 - 1 / (1 + exp(-k_global * (x_new[i] - x0_global))));
          y_hat_rel[i] = y_hat_new[i] / y_hat_new[1];

          y_prot_age[i, j] = 1 / (1 + exp(-k_temp[j] * (x_new[i] - x0_temp[j])));
          y_hat_new_age[i, j] =  L_temp[j] * (1 - 1 / (1 + exp(-k_temp[j] * (x_new[i] - x0_temp[j]))));
          y_hat_rel_age[i, j] = y_hat_new_age[i, j] / y_hat_new_age[1, j];
      }
    }
}