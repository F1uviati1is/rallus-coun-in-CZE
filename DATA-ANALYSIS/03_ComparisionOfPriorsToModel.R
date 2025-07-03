# Script for testing different model prior settings (UPDATED)

# SETUP ####
library(ggplot2)
library(glue)
library(rstan)
library(tibble)
library(tidyr)
library(dplyr)

# DATA ####
dta_dataForStanModel <- readRDS("DATA-ANALYSIS/data/res_dataStan.rds")

# FUNCTION TO GENERATE STAN CODE ####
fce_generateModelCode <- function(
    alpha_mean, alpha_sd, 
    beta_div_mean, beta_div_sd,
    beta_qual_mean = 0, beta_qual_sd,
    sigma_area_sd = 2,
    phi_meanlog, phi_sdlog
) {
  glue('
data {{
  int<lower=1> N;
  int<lower=1> J_area;
  int<lower=0> y[N];
  int<lower=1, upper=J_area> area[N];
  vector[N] environment_diversity;
  vector[N] environment_quality;
}}

parameters {{
  real<lower=0> alpha;
  real<lower=0> beta_diversity;
  real beta_quality;
  vector[J_area] u_area;
  real<lower=0> sigma_area;
  real<lower=0> phi;
}}

model {{
  alpha ~ lognormal({alpha_mean}, {alpha_sd});
  beta_diversity ~ lognormal({beta_div_mean}, {beta_div_sd});
  beta_quality ~ normal({beta_qual_mean}, {beta_qual_sd});
  u_area ~ normal(0, sigma_area);
  sigma_area ~ normal(0, {sigma_area_sd});
  phi ~ lognormal({phi_meanlog}, {phi_sdlog});

  for (n in 1:N) {{
    y[n] ~ neg_binomial_2_log(
      log(alpha) +
      beta_diversity * environment_diversity[n] +
      beta_quality * environment_quality[n] +
      u_area[area[n]],
      phi
    );
  }}
}}
')
}

# PRIOR SCENARIOS ####
dta_priorScenarios <- list(
  list(
    name = "infAlpha_weakOthers",
    alpha_mean = -1.3,
    alpha_sd = 0.5,
    beta_div_mean = 0,
    beta_div_sd = 4,
    beta_qual_sd = 4,
    phi_meanlog = log(1.62),
    phi_sdlog = 1
  ),
  list(
    name = "flatEverything",
    alpha_mean = -1,
    alpha_sd = 5,
    beta_div_mean = 0,
    beta_div_sd = 5,
    beta_qual_sd = 5,
    phi_meanlog = log(1.62),
    phi_sdlog = 1.5
  ),
  list(
    name = "tightAll",
    alpha_mean = -1.3,
    alpha_sd = 0.3,
    beta_div_mean = 0.5,
    beta_div_sd = 0.5,
    beta_qual_sd = 0.5,
    phi_meanlog = log(1.62),
    phi_sdlog = 0.2
  ),
  list(
    name = "shrinkCovariates",
    alpha_mean = -1.3,
    alpha_sd = 1,
    beta_div_mean = 0,
    beta_div_sd = 0.3,
    beta_qual_sd = 0.3,
    phi_meanlog = log(1.62),
    phi_sdlog = 0.5
  )
)

# RESULTS CONTAINER ####
fin_results <- list()

# MODEL FITTING ####
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
  
  set.seed(22)
  tmp_mod <- stan_model(model_code = tmp_code)
  
  set.seed(23)
  tmp_fit <- sampling(tmp_mod, data = dta_dataForStanModel, iter = 2000, chains = 4, seed = 23)
  
  fin_results[[scenario$name]] <- tmp_fit
}

# VISUALIZATION ####
alpha_samples_list <- lapply(names(fin_results), function(scenario_name) {
  alpha_samples <- rstan::extract(fin_results[[scenario_name]])$alpha
  tibble(
    scenario = scenario_name,
    alpha = alpha_samples * 11060
  )
})

dta_posteriorLong <- bind_rows(alpha_samples_list)

dta_means <- dta_posteriorLong |
  group_by(scenario) |
  summarise(mean_alpha = mean(alpha), .groups = "drop") |
  mutate(
    label = paste0("Mean = ", round(mean_alpha, 2)),
    y_position = 0.000004
  )

ggplot(dta_posteriorLong, aes(x = alpha, fill = scenario)) +
  geom_density(alpha = 0.5) +
  geom_vline(data = dta_means, aes(xintercept = mean_alpha, color = scenario),
             linetype = "dashed", linewidth = 1) +
  geom_text(
    data = dta_means,
    aes(x = mean_alpha, y = y_position, label = label),
    angle = 90,
    hjust = -0.1,
    size = 3.5,
    color = "black",     # <-- make text labels black
    show.legend = FALSE
  ) +
  labs(
    title = "Comparison of Posterior α × Area Across Priors",
    x = "Estimated Total Number of Individuals in CZ",
    y = "Density"
  ) +
  coord_cartesian(xlim = c(0, 25000), ylim = c(0, 3e-04)) +
  theme_minimal() +
  theme(
    text = element_text(color = "black"),           # Axis + title text
    axis.title = element_text(color = "black"),
    axis.text = element_text(color = "black"),
    legend.text = element_text(color = "black"),
    legend.title = element_text(color = "black"),
    plot.title = element_text(color = "black")
  )
