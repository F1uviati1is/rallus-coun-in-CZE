data {
  int<lower=1> N;                        // Number of observations
  int<lower=1> J_area;                   // Number of unique areas
  int<lower=0> y[N];                     // Response variable (count data)
  int<lower=1, upper=J_area> area[N];    // Area index
  vector[N] environment_diversity;       // Covariate: environmental diversity
  vector[N] environment_quality;         // Covariate: environmental quality
}

parameters {
  real<lower=0> alpha;                   // Baseline intercept (on original scale)
  real<lower=0> beta_diversity;          // Effect of diversity (positive constraint)
  real beta_quality;                     // Effect of quality
  vector[J_area] u_area;                 // Random effects for areas
  real<lower=0> sigma_area;              // SD of area random effects
  real<lower=0> phi;                     // Overdispersion parameter
}

model {
  // Flat (weakly informative) priors
  alpha ~ lognormal(-1, 5);
  beta_diversity ~ lognormal(0, 5);
  beta_quality ~ normal(0, 5);
  u_area ~ normal(0, sigma_area);
  sigma_area ~ normal(0, 2);
  phi ~ lognormal(log(1.62), 1.5);

  // Likelihood
  for (n in 1:N) {
    y[n] ~ neg_binomial_2_log(
      log(alpha) +
      beta_diversity * environment_diversity[n] +
      beta_quality * environment_quality[n] +
      u_area[area[n]],
      phi
    );
  }
}
