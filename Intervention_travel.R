# Load required libraries
# TRAVEL INTERVENTION
library(tidyverse)
library(brms)

# Ensure travel_frequency is numeric
data$travel_frequency <- as.numeric(as.character(data$travel_frequency))

# Set priors
priors <- get_custom_prior("bernoulli")

# 1. Fit observed model (total effect with covariates)
model_observed <- brm(
  infected_bin ~ travel_frequency + age_class + Sex_bin + Location,
  data = data,
  family = bernoulli(),
  prior = priors,
  iter = 10000, warmup = 3000, chains = 4, cores = 4,
  control = list(max_treedepth = 15)
)

# 2. Fit counterfactual model (total effect, backdoor paths are now closed so no confounders)
model_cf <- brm(
  infected_bin ~ travel_frequency,
  data = data,
  family = bernoulli(),
  prior = priors,
  iter = 10000, warmup = 3000, chains = 4, cores = 4,
  control = list(max_treedepth = 15)
)

# 3. Subset: individuals who currently travel daily (92 times)
daily_travellers <- data$travel_frequency == 92
data_daily <- data[daily_travellers, ] %>% droplevels()

# 4. Predict infection under observed model (original frequency)
obs_pred_daily <- rowMeans(posterior_epred(
  model_observed,
  newdata = data_daily,
  re_formula = NA
))

# 5. Create counterfactual: set travel_frequency = 0 for daily travellers
data_cf_daily <- data_daily %>%
  mutate(travel_frequency = 0) %>%
  droplevels()

# 6. Predict infection under counterfactual model
cf_pred_daily <- rowMeans(posterior_epred(
  model_cf,
  newdata = data_cf_daily,
  re_formula = NA
))

# 7. Summarise change in prevalence 
diff_daily <- cf_pred_daily - obs_pred_daily
quantile(diff_daily, probs = c(0.025, 0.5, 0.975))


# 8. Create counterfactual dataset for full population: set travel_frequency = 0
data_cf_all <- data %>%
  mutate(travel_frequency = 0) %>%
  droplevels()

# 9. Predict infection under observed model (original travel frequencies)
obs_pred_all <- rowMeans(posterior_epred(
  model_observed,
  newdata = data,
  re_formula = NA
))

# 10. Predict infection under counterfactual (everyone not travelling)
cf_pred_all <- rowMeans(posterior_epred(
  model_cf,
  newdata = data_cf_all,
  re_formula = NA
))

# 11. Summarise change in prevalence for the full population
diff_all <- cf_pred_all - obs_pred_all
quantile(diff_all, probs = c(0.025, 0.5, 0.975))

#====================================================================================
# Graphs     #######
# ===================================================================================
#
library(tidyverse)

#Daily traevllers only ####
#

# Prepare data
df_daily <- tibble(
  prob = c(obs_pred_daily, cf_pred_daily),
  type = rep(c("Observed", "Counterfactual"), each = length(obs_pred_daily))
)

# Plot
ggplot(df_daily, aes(x = prob, fill = type)) +
  geom_histogram(position = "identity", alpha = 0.6, bins = 30, colour = "black") +
  scale_fill_manual(values = c("Observed" = "#3B4BA3", "Counterfactual" = "#F4A261")) +
  labs(
    title = "Predicted Infection Probability (Daily Travellers)",
    x = "Infection probability",
    y = "Number of posterior predictions"
  ) +
  theme_minimal()

#Histogram for whole population

# Prepare data
df_all <- tibble(
  prob = c(obs_pred_all, cf_pred_all),
  type = rep(c("Observed", "Counterfactual"), each = length(obs_pred_all))
)

# Plot
ggplot(df_all, aes(x = prob, fill = type)) +
  geom_histogram(position = "identity", alpha = 0.6, bins = 30, colour = "black") +
  scale_fill_manual(values = c("Observed" = "#3B4BA3", "Counterfactual" = "#F4A261")) +
  labs(
    title = "Predicted Infection Probability (Whole Population)",
    x = "Infection probability",
    y = "Number of posterior predictions"
  ) +
  theme_minimal()+
  theme(text=element_text(size=14))

#side by side graphs
#
#
#library(tidyverse)
library(patchwork)  # for side-by-side plots

# Simulated draws (from your posterior_epred output)
# Make sure you have these ready:
# - obs_pred_all, cf_pred_all: full population (n_draws × n_individuals)
# - obs_pred_daily, cf_pred_daily: daily travellers (n_draws × n_individuals)

# Convert to long format
make_long_df <- function(pred_matrix, label) {
  as.data.frame(pred_matrix) %>%
    pivot_longer(everything(), names_to = "Individual", values_to = "Infection_Prob") %>%
    mutate(Type = label)
}

df_all_obs <- make_long_df(obs_pred_all, "Observed")
df_all_cf  <- make_long_df(cf_pred_all,  "Counterfactual")
df_daily_obs <- make_long_df(obs_pred_daily, "Observed")
df_daily_cf  <- make_long_df(cf_pred_daily,  "Counterfactual")

# Combine
df_all_combined   <- bind_rows(df_all_obs, df_all_cf)
df_daily_combined <- bind_rows(df_daily_obs, df_daily_cf)



# Plot settings
#use plot function from Intervention_activity.R called plot_hist2

p1 <- plot_hist2(df_all_combined,   "All individuals")
p2 <- plot_hist2(df_daily_combined, "Daily travellers")

# Side-by-side layout
side <-p1 + p2 + plot_layout(guides = "collect") & theme(legend.position = "bottom")

side

ggsave("Figures/Travel_intervention.png", side, width = 12, height = 10, dpi = 300)

#################
################## 