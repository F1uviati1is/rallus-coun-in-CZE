# Script for testing different model parameter settings

# SETUP ####
library(ggplot2)
library(glue)
library(rstan)
library(tibble)
library(tidyr)
library(dplyr)

# DATA ####
dta_dataForStanModel <- readRDS("DATA-ANALYSIS/data/res_dataStan.rds")

# FUNCTION: Generate model code with parameterized priors ####
fce_generateModelCode <- function(alpha_mean, alpha_sd, 
                                  beta_div_mean, beta_div_sd,
                                  beta_qual_sd,
                                  phi_meanlog, phi_sdlog) {
  glue('
data {{
  int<lower=1> N;
  int<lower=1> J_area;
  int<lower=1> T_time;
  int<lower=0> y[N];
  int<lower=1, upper=J_area> area[N];
  int<lower=1, upper=T_time> time[N];
  vector[N] environment_diversity;
  vector[N] environment_quality;
}}

parameters {{
  real<lower=0> alpha;
  real<lower=0> beta_diversity;
  real beta_quality;
  vector[J_area] u_area;
  vector[T_time] u_time;
  real<lower=0> sigma_area;
  real<lower=0> sigma_time;
  real<lower=0> phi;
}}

model {{
  alpha ~ lognormal({alpha_mean}, {alpha_sd});
  beta_diversity ~ lognormal({beta_div_mean}, {beta_div_sd});
  beta_quality ~ normal(0, {beta_qual_sd});
  u_area ~ normal(0, sigma_area);
  u_time ~ normal(0, sigma_time);
  sigma_area ~ normal(0, 2);
  sigma_time ~ normal(0, 2);
  phi ~ lognormal({phi_meanlog}, {phi_sdlog});

  for (n in 1:N) {{
    y[n] ~ neg_binomial_2_log(
      log(alpha) +
      beta_diversity * environment_diversity[n] +
      beta_quality * environment_quality[n] +
      u_area[area[n]] +
      u_time[time[n]],
      phi
    );
  }}
}}
')
}

# LIST OF PRIOR SCENARIOS ####
dta_priorScenarios <- list(
  list(name = "default",
       alpha_mean = -5.986,
       alpha_sd = 3,
       beta_div_mean = 0,
       beta_div_sd = 2,
       beta_qual_sd = 2,
       phi_meanlog = log(1.62),
       phi_sdlog = 0.5),
  
  list(name = "weaker_alpha",
       alpha_mean = -3,
       alpha_sd = 1,
       beta_div_mean = 0,
       beta_div_sd = 2,
       beta_qual_sd = 2,
       phi_meanlog = log(1.62),
       phi_sdlog = 0.5)
)

# STORAGE FOR RESULTS ####
fin_results <- list()

# MODELLING ####
for (scenario in dta_priorScenarios) {
  cat("Running:", scenario$name, "\n")
  
  tmp_code <- fce_generateModelCode(
    alpha_mean = scenario$alpha_mean,
    alpha_sd = scenario$alpha_sd,
    beta_div_mean = scenario$beta_div_mean,
    beta_div_sd = scenario$beta_div_sd,
    beta_qual_sd = scenario$beta_qual_sd,
    phi_meanlog = scenario$phi_meanlog,
    phi_sdlog = scenario$phi_sdlog
  )
  
  tmp_mod <- stan_model(model_code = tmp_code)
  
  set.seed(23)
  tmp_fit <- sampling(tmp_mod, data = dta_dataForStanModel, iter = 2000, chains = 4, seed = 23)
  
  fin_results[[scenario$name]] <- tmp_fit
}

# VISUALIZATION ####
# Extract and combine alpha samples from each scenario
alpha_samples_list <- lapply(names(fin_results), function(scenario_name) {
  alpha_samples <- rstan::extract(fin_results[[scenario_name]])$alpha
  tibble(
    scenario = scenario_name,
    alpha = alpha_samples * 11060
  )
})

dta_posteriorLong <- bind_rows(alpha_samples_list)

# Compute means for each scenario
dta_means <- dta_posteriorLong |>
  group_by(scenario) |>
  summarise(mean_alpha = mean(alpha), .groups = "drop") |>
  mutate(
    label = paste0("Mean = ", round(mean_alpha, 2)),
    y_position = 0.000004  # adjust this depending on your density height
  )

# Plot
ggplot(dta_posteriorLong, aes(x = alpha, fill = scenario)) +
  geom_density(alpha = 0.5) +
  geom_vline(data = dta_means, aes(xintercept = mean_alpha, color = scenario),
             linetype = "dashed", linewidth = 1) +
  geom_text(
    data = dta_means,
    aes(x = mean_alpha, y = y_position, label = label, color = scenario),
    angle = 90,
    hjust = -0.1,
    size = 3.5,
    show.legend = FALSE
  ) +
  labs(
    title = "Comparison of Posterior α × Area Across Priors",
    x = "Estimated Total Number of Individuals in CZ",
    y = "Density"
  ) +
  theme_minimal()
