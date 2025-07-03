## SETUP ####

library(bayesplot)
library(dplyr)
library(ggplot2)
library(MASS)
library(rstan)
library(purrr)

## DATA LOAD ####

# Load data list
dta_dataForStanModel <- readRDS("DATA-ANALYSIS/data/res_dataStan.rds")

# Load Stan models
mod_model0 <- readLines("DATA-ANALYSIS/data/res_model0.stan")
mod_model1 <- readLines("DATA-ANALYSIS/data/res_model1.stan")

### COMPILING THE MODELS ####

stan_model0 <- stan_model(model_code = mod_model0)
stan_model1 <- stan_model(model_code = mod_model1)


set.seed(23)
fit0 <- sampling(
  stan_model0,
  data = dta_dataForStanModel,
  iter = 2000,
  chains = 4
)

set.seed(23)
fit1 <- sampling(
  stan_model1,
  data = dta_dataForStanModel,
  iter = 2000,
  chains = 4
)


#### VISUALIZATION ####
# Extract alpha samples for both models (and multiply by proportional mapped area to whole CZE)
alpha_samples0 <- extract(fit0)$alpha * 11060
alpha_samples1 <- extract(fit1)$alpha * 11060

# Create a combined dataframe
alpha_samples_df <- data.frame(
  samples = c(alpha_samples0, alpha_samples1),
  model = rep(c("Model 0", "Model 1"), each = length(alpha_samples0))
)

# Compute statistics for annotations
mean_alpha0 <- mean(alpha_samples0)
q1_alpha0 <- quantile(alpha_samples0, 0.25)
q3_alpha0 <- quantile(alpha_samples0, 0.75)

mean_alpha1 <- mean(alpha_samples1)
q1_alpha1 <- quantile(alpha_samples1, 0.25)
q3_alpha1 <- quantile(alpha_samples1, 0.75)

# Plot both distributions with limited x-axis
ggplot(alpha_samples_df, aes(x = samples, fill = model)) +
  geom_density(alpha = 0.3, show.legend = TRUE) +  # Enable legend for distinction
  scale_fill_manual(values = c("yellow", "green")) +  # Model 0 in yellow, Model 1 in green
  geom_vline(aes(xintercept = mean_alpha0), color = "yellow", linetype = "dashed", linewidth = 1) +
  geom_text(aes(x = mean_alpha0, y = 0, label = paste("Mean model 0 =", round(mean_alpha0, 2))), color = "blue", vjust = -1) +
  geom_vline(aes(xintercept = mean_alpha1), color = "green", linetype = "dashed", linewidth = 1) +
  geom_text(aes(x = mean_alpha1, y = 0, label = paste("Mean model 1 =", round(mean_alpha1, 2))), color = "blue", vjust = -3) +
  labs(title = "",
       x = "Počet chřástalů",
       y = "Pravděpodobnost") +
  theme_minimal() +
  xlim(0, 20000)  # Limit x-axis to 0 - 20000


#### LOLO-CV MODEL DIAGNOSTICS #### (Leave one locality out - cross validation)


