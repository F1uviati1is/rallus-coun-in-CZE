# For preparation of data

## SETUP ####

library(readxl)
library(sf)
library(dplyr)
library(vegan)

## DATA LOAD ####
# Rail data per point
org_railPerPoint <- st_read("DATA-ANALYSIS/data/body_chrastal.shp", options = "ENCODING=UTF-8")
dta_railCounts <- org_railPerPoint

# Naming columns to chose data from to count Shanon index
tmp_cols <- c("rakos", "orobinec", "zblochan", "ostrice", "dreviny", "ostatni")

## DATA EDITING ####
# This following awefully long data managing is, selecting right methodology, cutting useless columns, computing Shanon index, creating Termin variable, computing usage of number of rails for the likelihood
dta_railCounts <- dta_railCounts |>
  dplyr::filter(Metoda %in% c("bodový transekt - intenzivní", "AKU - intenzivní")) |>  # Counting only specific methods
  dplyr::select(Name, Oblast, X1_termin, X2_termin, X3_termin, zarust, ostatni, rakos, orobinec, zblochan, ostrice, dreviny, Celkem)|> # Using only needed columns
  mutate(
    Shanon = rowSums(across(all_of(tmp_cols), ~ .x / (.x + 1))), # Shanon index
    Shanon_vegan = diversity(as.matrix(across(all_of(tmp_cols))), index = "shannon"), # Shanon z veganu při 100% hodnotě jediného biotopu vkládá 0, zatímco sloupec Shanon (vzorec který nepřevádí čísla na přirozené logaritmy než je dělí) vkládá hodnotu as 0.9 - To mi přijde lepší, protože pro celou lokalitu to taky znamená nějakou diverzitu pokud je tam ač jediný biotop.
    pocet_biotopu = rowSums(across(all_of(tmp_cols), ~ .x > 0)) # Number of biotopes
  )|>
  group_by(Oblast) |>
  mutate(
    sum_X1 = sum(X1_termin, na.rm = TRUE),
    sum_X2 = sum(X2_termin, na.rm = TRUE),
    sum_X3 = sum(X3_termin, na.rm = TRUE),
    Termin = case_when(                          # Creating Termin Column with number of termin used for model
      sum_X1 > sum_X2 & sum_X1 > sum_X3 ~ 1,
      sum_X2 > sum_X1 & sum_X2 > sum_X3 ~ 2,
      sum_X3 > sum_X1 & sum_X3 > sum_X2 ~ 3,
      TRUE ~ 1
    ),
    water_rails_count = case_when(               # Creating alpha! number of rails taken from a Termin of most rails counted
      Termin == 1 ~ X1_termin,
      Termin == 2 ~ X2_termin,
      Termin == 3 ~ X3_termin,
      TRUE ~ X1_termin
    )
  ) |>
  ungroup() |>
  mutate(
    water_rails_count = case_when(                # For selected Oblast use Celkem column of number of rails, because in data it is not written in Termines (coz only one counting)
      Oblast %in% c("Třebíčsko", "Telčsko", "Prostějovsko") ~ Celkem,
      TRUE ~ water_rails_count
    )
  )

dta_railCounts$zarust <- as.numeric(scale(dta_railCounts$zarust))   # Scaling explanatory variable coz it was making mess
dta_railCounts$Termin <- as.numeric(scale(dta_railCounts$Termin))   # Scaling explanatory variable coz it was making mess

# Stan data
fin_dataStan <- list(
  N = nrow(dta_railCounts),
  J_area = length(unique(dta_railCounts$Oblast)),
  T_time = length(unique(dta_railCounts$Termin)),  
  y = dta_railCounts$water_rails_count,
  area = as.numeric(factor(dta_railCounts$Oblast)),
  time = as.numeric(factor(dta_railCounts$Termin)),
  environment_diversity = dta_railCounts$Shanon,
  environment_quality = dta_railCounts$zarust  
)

## MODEL WRITING ####

# Model
fin_model0 <- "data {
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
  alpha ~ lognormal(-5.986, 3);      // Prior for baseline
  u_area ~ normal(0, sigma_area);    // Prior for area-specific variation
  u_time ~ normal(0, sigma_time);    // Prior for time-specific variation
  sigma_area ~ normal(0, 2);         // Prior for area variation
  sigma_time ~ normal(0, 2);         // Prior for time variation
  phi ~ lognormal(log(1.62), 0.5);   // Prior for overdispersion
  
  // Likelihood
  for (n in 1:N) {
    y[n] ~ neg_binomial_2_log(
      log(alpha) +
        u_area[area[n]] +  // Random effect for area
      u_time[time[n]],   // Random effect for time
      phi                // Overdispersion parameter
    );
  }
}"

fin_model_text1 <- "
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
