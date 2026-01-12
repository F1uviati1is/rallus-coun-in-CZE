data {
  int<lower=1> N;                   // Number of observations
  int<lower=1> J_area;              // Number of unique areas
  int<lower=1> T_time;              // Number of unique time points
  int<lower=0> y[N];                // Response variable (count data)
  int<lower=1, upper=J_area> area[N]; // Area index for each observation
  int<lower=1, upper=T_time> time[N]; // Time index for each observation
}

parameters {
  real<lower=0> alpha;              // Baseline parameter (intercept)
  vector[J_area] u_area;            // Random effects for areas
  vector[T_time] u_time;            // Random effects for time points
  real<lower=0> sigma_area;         // Standard deviation of area effects
  real<lower=0> sigma_time;         // Standard deviation of time effects
  real<lower=0> phi;                // Overdispersion parameter (negative binomial)
}

model {
  // Priors
  alpha ~ lognormal(-1, 5);          // Prior for baseline
  beta_diversity ~ lognormal(0,5);
  beta_quality ~ normal(0,5);
  u_area ~ normal(0, sigma_area);    // Prior for area-specific variation
  u_time ~ normal(0, sigma_time);    // Prior for time-specific variation
  sigma_area ~ normal(0, 2);         // Prior for area variation
  sigma_time ~ normal(0, 2);         // Prior for time variation
  phi ~ lognormal(log(1.62), 1.5);   // Prior for overdispersion
  
  // Likelihood
  for (n in 1:N) {
    y[n] ~ neg_binomial_2_log(
      log(alpha) +
        u_area[area[n]] +  // Random effect for area
      u_time[time[n]],   // Random effect for time
      phi                // Overdispersion parameter
    );
  }
}
