# Script of testing different model aprameters setting

# SETUP ####

library(ggplot2)
library(rstan)

# DATA ####
dta_dataForStanModel <- readRDS("DATA-ANALYSIS/data/res_dataStan.rds")

mod_test <- "
data {
  int<lower=1> N;                   // Number of observations
  int<lower=1> J_area;              // Number of unique areas
  int<lower=1> T_time;              // Number of unique time points
  int<lower=0> y[N];                // Response variable (count data)
  int<lower=1, upper=J_area> area[N]; // Area index for each observation
  int<lower=1, upper=T_time> time[N]; // Time index for each observation
  vector[N] environment_diversity;  // Predictor: environmental diversity
  vector[N] environment_quality;    // Predictor: environmental quality
}

parameters {
  real<lower=0> alpha;              // Baseline parameter (intercept)
  real<lower=0> beta_diversity;     // Effect of environmental diversity
  real beta_quality;                // Effect of environmental quality
  vector[J_area] u_area;            // Random effects for areas
  vector[T_time] u_time;            // Random effects for time points
  real<lower=0> sigma_area;         // Standard deviation of area effects
  real<lower=0> sigma_time;         // Standard deviation of time effects
  real<lower=0> phi;                // Overdispersion parameter (negative binomial)
}

model {
  // Priors
  alpha ~ lognormal(-5.986, 3);      // Prior for baseline (0.226 * 11060 = 2499.56)
  beta_diversity ~ lognormal(0,2);   // Prior for diversity effect
  beta_quality ~ normal(0, 2);       // Prior for quality effect
  u_area ~ normal(0, sigma_area);    // Prior for area-specific variation
  u_time ~ normal(0, sigma_time); // Prior for time-specific variation
  sigma_area ~ normal(0, 2);         // Normal prior for area variation
  sigma_time ~ normal(0, 2);         // Normal prior for time variation
  phi ~ lognormal(log(1.62), 0.5);;     // Truncated normal prior for overdispersion

  // Likelihood
  for (n in 1:N) {
    y[n] ~ neg_binomial_2_log( // Negative binomial model with log-link
      log(alpha) +             // Baseline count (in log space)
      beta_diversity * environment_diversity[n] + // Effect of diversity
      beta_quality * environment_quality[n] +     // Effect of quality
      u_area[area[n]] +  // Random effect for area
      u_time[time[n]],   // Random effect for time
      phi                // Overdispersion parameter
    );
  }
}
"

# MODEL

mod_test <- stan_model(model_code = mod_test)

set.seed(23)
fit_modTest <- sampling(
  mod_test,
  data = dta_dataForStanModel,
  iter = 2000,
  chains = 4
)

#### VISUALIZATION ####
# Extract alpha samples for both models (and multiply by proportional mapped area to whole CZE)
dta_alphaSamples <- extract(fit_modTest)$alpha * 11060

# Compute statistics for annotations
tmp_meanAlpha <- mean(dta_alphaSamples)
tmp_q1Alpha <- quantile(dta_alphaSamples, 0.25)
tmp_q3Alpha <- quantile(dta_alphaSamples, 0.75)

# Convert to data frame for ggplot
dta_alphaSamples_df <- data.frame(alpha = dta_alphaSamples)

# Plot
ggplot(dta_alphaSamples_df, aes(x = alpha)) +
  geom_density(fill = "steelblue", alpha = 0.6) +
  geom_vline(xintercept = tmp_meanAlpha, color = "red", linetype = "dashed", size = 1) +
  geom_vline(xintercept = tmp_q1Alpha, color = "darkgreen", linetype = "dotted", size = 1) +
  geom_vline(xintercept = tmp_q3Alpha, color = "darkgreen", linetype = "dotted", size = 1) +
  labs(
    title = "Posterior Distribution of Total Individuals (α × Area)",
    x = "Estimated Total Number of Individuals in CZ",
    y = "Density"
  ) +
  theme_minimal()
