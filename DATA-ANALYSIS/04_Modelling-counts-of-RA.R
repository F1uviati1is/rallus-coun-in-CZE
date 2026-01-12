# For final counting the model of rails, visualize and diagnose

## SETUP ####

library(bayesplot)
library(dplyr)
library(ggplot2)
library(MASS)
library(rstan)
library(purrr)
rm(list = ls())

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
  chains = 4,
  control = list(adapt_delta = 0.95)  #Coz 
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

#### DIAGNOSTICSS ####

# Table
rstan::monitor(fit, warmup = 1000, print = TRUE)

# Traceplot
traceplot(fit, pars = c("alpha", "beta_diversity", "beta_quality", "phi", "sigma_area"))

check_hmc_diagnostics(fit)

post <- rstan::extract(fit)
y_rep <- replicate(100, {
  rnbinom(length(dta_dataForStanModel$y),
          mu = exp(log(post$alpha[1]) + 
                     post$beta_diversity[1] * dta_dataForStanModel$environment_diversity +
                     post$beta_quality[1] * dta_dataForStanModel$environment_quality +
                     post$u_area[1, dta_dataForStanModel$area]),
          size = post$phi[1])
})

y_rep <- t(y_rep)  # Transpose to [100 draws × 885 observations]

bayesplot::ppc_dens_overlay(dta_dataForStanModel$y, y_rep)

posterior <- as.data.frame(fit)
ggplot(posterior, aes(x = beta_diversity)) + geom_density(fill = "skyblue")
ggplot(posterior, aes(x = beta_quality)) + geom_density(fill = "skyblue")
ggplot(posterior, aes(x = phi)) + geom_density(fill = "skyblue")


#### LOLO-CV MODEL DIAGNOSTICS #### (Leave one locality out - cross validation)


