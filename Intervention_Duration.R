# Load required libraries
library(tidyverse)
library(brms)

#=================================================================
#DURATION INTERVENTION########################
#================================================================
# Set priors
priors <- get_custom_prior("bernoulli")



# 1. Fit observed model (total effect with covariates)
model_observed_dur <- brm(
  infected_bin ~ DurMin +age_class+Location+Sex_bin+ActNEW,
  data = data,
  family = bernoulli(),
  prior = priors,
  iter = 10000, warmup = 3000, chains = 4, cores = 4,
  control = list(max_treedepth = 15)
)

# 2. Fit counterfactual model (total effect, backdoor paths are now closed so no confounders)
model_cf_dur <- brm(
  infected_bin ~ DurMin,
  data = data,
  family = bernoulli(),
  prior = priors,
  iter = 10000, warmup = 3000, chains = 4, cores = 4,
  control = list(max_treedepth = 15)
)


# 3. Predict infection under observed model (original frequency)
obs_pred_dur <- rowMeans(posterior_epred(
  model_observed_dur,
  newdata = data,
  re_formula = NA
))

# 4. Create counterfactual: no duration in the lake for everyone
data_cf_du0<- data %>%
  mutate(DurMin = case_when(DurMin >=0 ~0)) %>%
  droplevels()

# 5. Predict infection under counterfactual model
cf_pred_dur <- rowMeans(posterior_epred(
  model_cf_dur,
  newdata = data_cf_du0,
  re_formula = NA
))

# 6. Summarise change in prevalence 
diff_dur<- cf_pred_dur - obs_pred_dur
quantile(diff_dur, probs = c(0.025, 0.5, 0.975))

#Now do it with only people who reported spending any time in the lake


data2 <- data %>%
  filter(DurMin >0)

# 1. Fit observed model (total effect with covariates)
model_observed_dur2 <- brm(
  infected_bin ~ DurMin +age_class+Location+Sex_bin+ActNEW,
  data = data2,
  family = bernoulli(),
  prior = priors,
  iter = 10000, warmup = 3000, chains = 4, cores = 4,
  control = list(max_treedepth = 15)
)

# 2. Fit counterfactual model (total effect, backdoor paths are now closed so no confounders)
model_cf_dur2 <- brm(
  infected_bin ~ DurMin,
  data = data2,
  family = bernoulli(),
  prior = priors,
  iter = 10000, warmup = 3000, chains = 4, cores = 4,
  control = list(max_treedepth = 15)
)


# 3. Predict infection under observed model (original frequency)
obs_pred_dur2 <- rowMeans(posterior_epred(
  model_observed_dur2,
  newdata = data2,
  re_formula = NA
))

# 4. Create counterfactual: no duration in the lake for everyone
data_cf_du02<- data2 %>%
  mutate(DurMin = case_when(DurMin >=0 ~0)) %>%
  droplevels()

# 5. Predict infection under counterfactual model
cf_pred_dur2 <- rowMeans(posterior_epred(
  model_cf_dur2,
  newdata = data_cf_du02,
  re_formula = NA
))

# 6. Summarise change in prevalence 
diff_dur2<- cf_pred_dur2 - obs_pred_dur2
quantile(diff_dur2, probs = c(0.025, 0.5, 0.975))


#====================================================================================
# Graphs     #######
# ===================================================================================
#
library(tidyverse)

 ####
#Everyone in the population

# Prepare data
df_dur_ALL <- tibble(
  prob = c(obs_pred_dur, cf_pred_dur),
  type = rep(c("Observed", "Counterfactual"), each = length(obs_pred_dur))
)

# Plot
ggplot(df_dur_ALL, aes(x = prob, fill = type)) +
  geom_histogram(position = "identity", alpha = 0.6, bins = 30, colour = "black") +
  scale_fill_manual(values = c("Observed" = "#3B4BA3", "Counterfactual" = "#F4A261")) +
  labs(
    title = "Predicted Infection Probability (Whole sample population)",
    x = "Infection probability",
    y = "Number of posterior predictions"
  ) +
  theme_minimal()

#Histogram foronly people who report contact with the lake

# Prepare data
df_dur <- tibble(
  prob = c(obs_pred_dur2, cf_pred_dur2),
  type = rep(c("Observed", "Counterfactual"), each = length(obs_pred_dur2))
)

# Plot
ggplot(df_dur, aes(x = prob, fill = type)) +
  geom_histogram(position = "identity", alpha = 0.6, bins = 30, colour = "black") +
  scale_fill_manual(values = c("Observed" = "#3B4BA3", "Counterfactual" = "#F4A261")) +
  labs(
    title = "Predicted Infection Probability (Lake contacters)",
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
make_long_df_dur <- function(pred_matrix, label) {
  as.data.frame(pred_matrix) %>%
    pivot_longer(everything(), names_to = "Individual", values_to = "Infection_Prob") %>%
    mutate(Type = label)
}

df_all_obs_dur <- make_long_df(obs_pred_dur, "Observed")
df_all_cf_dur  <- make_long_df(cf_pred_dur,  "Counterfactual")
df_dur2_obs <- make_long_df(obs_pred_dur2, "Observed")
df_dur2_cf  <- make_long_df(cf_pred_dur2,  "Counterfactual")

# Combine
df_all_combined_dur   <- bind_rows(df_all_obs_dur, df_all_cf_dur)
df_dur2_combined <- bind_rows(df_dur2_obs, df_dur2_cf)

# Plot settings
#use plot function from Intervention_activity.R called plot_hist2

p1 <- plot_hist2(df_all_combined_dur,   "All individuals")
p2 <- plot_hist2(df_dur2_combined, "Lake contactors")

# Side-by-side layout
side_dur <-p1 + p2 + plot_layout(guides = "collect") & theme(legend.position = "bottom")

side_dur

ggsave("Figures/Duration_intervention.png", side_dur, width = 12, height = 10, dpi = 300)

#################
################## 
#What about on worm burden

