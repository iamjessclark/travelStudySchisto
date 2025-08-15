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
data$travel_frequency <- as.factor(data$travel_frequency)

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

# Flag daily travellers in the factual data (before any CF edits)
data <- data %>%
  mutate(was_daily = travel_frequency == "92")

# Save factor levels so we can reapply them after mutate
tf_levels <- levels(data$travel_frequency)
act_levels <- if (is.factor(data$ActNEW)) levels(data$ActNEW) else NULL

cf_data_daily <- data %>%
  mutate(
    travel_frequency = if_else(travel_frequency == "92", "0", as.character(travel_frequency)),
    DurMin = if_else(travel_frequency == "0", 0, DurMin),
    ActNEW = if_else(travel_frequency == "0", "None", as.character(ActNEW))
  ) %>%
  mutate(
    travel_frequency = factor(travel_frequency, levels = tf_levels),
    ActNEW = if (!is.null(act_levels)) factor(ActNEW, levels = act_levels) else ActNEW
  )



# Population-wide infection probability per draw
# factual
obs_all_draws <- (posterior_epred(model_observed, 
                                         newdata = data, 
                                         re_formula = NA))
# counterfactual
cf_all_draws  <- (posterior_epred(model_observed, 
                                         newdata = cf_data_daily, 
                                         re_formula = NA))




#populstion wide mean per draw and CrIs
#
obs_mean_all <- rowMeans(obs_all_draws)
cf_mean_all  <- rowMeans(cf_all_draws)

# DAILY TRAVELLERS ONLY (use the flag to pick columns)
was_daily <- data$was_daily
stopifnot(any(was_daily))  # sanity check

obs_daily_draws <- rowMeans(obs_all_draws[, was_daily, drop = FALSE], na.rm = TRUE)
cf_daily_draws  <- rowMeans(cf_all_draws[,  was_daily, drop = FALSE], na.rm = TRUE)


#a function to get the median and crIs
summ <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]                 # drop NA / NaN / Inf
  c(
    median = stats::median(x),
    l95    = stats::quantile(x, 0.025, names = FALSE),
    u95    = stats::quantile(x, 0.975, names = FALSE)
  )
}


pop_obs_summary  <- 100 * summ(obs_mean_all)
pop_cf_summary   <- 100 * summ(cf_mean_all)
pop_diff_summary <- 100 * summ(cf_mean_all - obs_mean_all)



daily_obs_summary  <- 100 * summ(obs_daily_draws)
daily_cf_summary   <- 100 * summ(cf_daily_draws)
daily_diff_summary <- 100 * summ(cf_daily_draws - obs_daily_draws)

#put in table I can quote in text

library(tibble)
results_tbl <- tibble(
  group      = rep(c("Whole population", "Daily travellers"), each = 3),
  scenario   = rep(c("Observed", "Counterfactual", "Difference (CF−Obs)"), times = 2),
  median_pct = c(pop_obs_summary["median"],  pop_cf_summary["median"],  pop_diff_summary["median"],
                 daily_obs_summary["median"], daily_cf_summary["median"], daily_diff_summary["median"]),
  l95_pct    = c(pop_obs_summary["l95"],     pop_cf_summary["l95"],     pop_diff_summary["l95"],
                 daily_obs_summary["l95"],    daily_cf_summary["l95"],    daily_diff_summary["l95"]),
  u95_pct    = c(pop_obs_summary["u95"],     pop_cf_summary["u95"],     pop_diff_summary["u95"],
                 daily_obs_summary["u95"],    daily_cf_summary["u95"],    daily_diff_summary["u95"])
)
print(results_tbl)






str(data$travel_frequency); str(cf_data_daily$travel_frequency)
setdiff(levels(cf_data_daily$travel_frequency), levels(data$travel_frequency))
sum(is.na(obs_all_draws)); sum(is.na(cf_all_draws))



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
  prob = c(obs_daily_draws, cf_daily_draws),
  type = rep(c("Observed", "Counterfactual"), each = length(obs_daily_draws))
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


