# For final counting the model of rails, visualize and diagnose

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
mod_railPerPoint <- readLines("DATA-ANALYSIS/data/res_model.stan")

### COMPILING THE MODELS ####

stan_model <- stan_model(model_code = mod_railPerPoint)

set.seed(23)
fit <- sampling(
  stan_model,
  data = dta_dataForStanModel,
  iter = 2000,
  chains = 4
)

#### VISUALIZATION ####
# Extract alpha samples for both models (and multiply by proportional mapped area to whole CZE)
alpha_samples <- rstan::extract(fit)$alpha * 11060

# Create a combined dataframe
alpha_samples_df <- data.frame(
  samples = alpha_samples)

# Compute statistics for annotations
mean_alpha <- mean(alpha_samples)
q1_alpha <- quantile(alpha_samples, 0.25)
q3_alpha <- quantile(alpha_samples, 0.75)

# Plot both distributions with limited x-axis
ggplot(alpha_samples_df, aes(x = samples, fill = T)) +
  geom_density(alpha = 0.3, show.legend = F) +  # Enable legend for distinction
  scale_fill_manual(values = "yellow") +  
  geom_vline(aes(xintercept = mean_alpha), color = "black", linetype = "dashed", linewidth = 1) +
  geom_text(aes(x = mean_alpha, y = 0, label = paste("Mean =", round(mean_alpha, 2))), color = "blue", vjust = -1) +
  labs(title = "",
       x = "Počet chřástalů",
       y = "Pravděpodobnost") +
  theme_minimal() +
  xlim(0, 20000) # Limit x-axis to 0 - 20000


#### LOLO-CV MODEL DIAGNOSTICS #### (Leave one locality out - cross validation)


