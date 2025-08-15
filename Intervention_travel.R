# Load required libraries
# TRAVEL INTERVENTION
rm(list = ls())
library(tidyverse)
library(brms)



#install data
#
data <-read.csv("Travel_prevEDIT.csv")
nrow(data) #586

# Ensure travel_frequency is numeric
data$travel_frequency <- as.numeric(as.character(data$travel_frequency))

# Set priors
priors <- get_custom_prior("bernoulli")

# 1. Fit observed model
model_observed <- brm(
  infected_bin ~ travel_frequency + age_class + Sex_bin + Location,
  data = data,
  family = bernoulli(),
  prior = priors,
  iter = 10000, warmup = 3000, chains = 4, cores = 4,
  control = list(max_treedepth = 15)
)

#2. Set up counterfactual data
# daily travellers have zero travel and the downstream effects of this are also noted
# anyone who hase treavl = zero must also have duration = 0 and activity = none
cf_data_daily <- data %>%
  mutate(
    travel_frequency = if_else(travel_frequency == 92, 0, travel_frequency),
    DurMin = if_else(travel_frequency == 0, 0, DurMin),
    ActNEW = if_else(travel_frequency == 0, "None", ActNEW)
  )


# Population-wide mean infection probability per draw
obs_mean_all <- rowMeans(posterior_epred(model_observed, newdata = data, re_formula = NA))
cf_mean_all  <- rowMeans(posterior_epred(model_observed, newdata = cf_data_daily, re_formula = NA))

# Difference in percentage points
diff_all <- cf_mean_all - obs_mean_all
quantile(diff_all, probs = c(0.025, 0.5, 0.975))

# ---- Daily travellers only ----
daily_index <- which(data$travel_frequency == 92)

# Means within daily travellers group
obs_mean_daily <- rowMeans(posterior_epred(model_observed, newdata = data[daily_index, ], re_formula = NA))
cf_mean_daily  <- rowMeans(posterior_epred(model_observed, newdata = cf_data_daily[daily_index, ], re_formula = NA))

# Difference for daily travellers only
diff_daily <- cf_mean_daily - obs_mean_daily
quantile(diff_daily, probs = c(0.025, 0.5, 0.975))

#make two graphs

library(ggplot2)
library(tibble)
library(dplyr)

# --- Histogram for population-wide effect ---
df_all <- tibble(
  prob = c(obs_mean_all, cf_mean_all),
  type = rep(c("Observed", "Counterfactual"), each = length(obs_mean_all))
)

ggplot(df_all, aes(x = prob, fill = type)) +
  geom_histogram(position = "identity", alpha = 0.6, bins = 30, colour = "black") +
  scale_fill_manual(values = c("Observed" = "#3B4BA3", "Counterfactual" = "#F4A261")) +
  labs(
    title = "Predicted Infection Probability (Whole Population)",
    x = "Mean infection probability",
    y = "Number of posterior draws"
  ) +
  theme_minimal(base_size = 14)

# --- Histogram for daily travellers only ---
df_daily <- tibble(
  prob = c(obs_mean_daily, cf_mean_daily),
  type = rep(c("Observed", "Counterfactual"), each = length(obs_mean_daily))
)

ggplot(df_daily, aes(x = prob, fill = type)) +
  geom_histogram(position = "identity", alpha = 0.6, bins = 30, colour = "black") +
  scale_fill_manual(values = c("Observed" = "#3B4BA3", "Counterfactual" = "#F4A261")) +
  labs(
    title = "Predicted Infection Probability (Daily Travellers)",
    x = "Mean infection probability",
    y = "Number of posterior draws"
  ) +
  theme_minimal(base_size = 14)


