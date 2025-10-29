
rm(list = ls())

# Load required libraries
library(tidyverse)
library(brms)

#=================================================================
#DURATION INTERVENTION########################
#================================================================
# Set priors
priors <- get_custom_prior("bernoulli")

#install data
#
data <-read.csv("Travel_prevEDIT.csv")
nrow(data) #586



# 1. Fit observed model (total effect with covariates)
model_observed_dur <- brm(
  infected_bin ~ DurMin +age_class+Location+Sex_bin+ActNEW,
  data = data,
  family = bernoulli(),
  prior = priors,
  iter = 10000, warmup = 3000, chains = 4, cores = 4,
  control = list(max_treedepth = 15)
)

#2. make counterfactual dataset where time in water is set to 0

cf_data_dur <-  data %>%
  mutate(DurMin = 0)

cf_data_dur
  
# 3. Predict infection under observed model (original frequency)
obs_pred_dur <- posterior_epred(
  model_observed_dur,
  newdata = data,
  re_formula = NA
)


# 4 Predict infection under counterfactual model
cf_pred_dur <- posterior_epred(
  model_observed_dur,
  newdata = cf_data_dur,
  re_formula = NA
)

############
############
# --- after you have obs_pred_dur and cf_pred_dur ---

data$travel_frequency <- as.factor(data$travel_frequency)
# helper to summarise a vector of draws
summ <- function(x) c(median = median(x),
                      l95 = quantile(x, 0.025),
                      u95 = quantile(x, 0.975))

# function: per-draw group means (pred_matrix: draws x individuals)
# robust: always returns a numeric matrix (draws x levels)
group_means <- function(pred_matrix, groups_factor) {
  levs <- levels(groups_factor)
  n_draws <- nrow(pred_matrix)
  
  out <- vapply(
    levs,
    function(lv) {
      idx <- which(groups_factor == lv)
      if (length(idx) == 0) rep(NA_real_, n_draws)
      else rowMeans(pred_matrix[, idx, drop = FALSE])
    },
    numeric(n_draws)  # enforce numeric
  )
  colnames(out) <- levs
  out
}

obs_pred_dur <- as.matrix(obs_pred_dur)
storage.mode(obs_pred_dur) <- "double"

cf_pred_dur  <- as.matrix(cf_pred_dur)
storage.mode(cf_pred_dur)  <- "double"

# per-draw mean risk within each travel_frequency level
obs_by_freq <- group_means(obs_pred_dur, data$travel_frequency)
cf_by_freq  <- group_means(cf_pred_dur,  data$travel_frequency)

# per-draw contrasts (CF − OBS) within each level = CATEs
diff_by_freq <- cf_by_freq - obs_by_freq

# summarise for each level (risk difference on probability scale)
cate_table <- t(apply(diff_by_freq, 2, function(col) summ(col)))
cate_table

############ mechaistic counterfactual

#populstion wide mean per draw and CrIs
#
obs_mean_all <- rowMeans(obs_pred_dur)
cf_mean_all  <- rowMeans(cf_pred_dur)




# 6. Summarise change in prevalence 
diff_dur<- cf_mean_all - obs_mean_all
quantile(diff_dur, probs = c(0.025, 0.5, 0.975))




summ <- function(x) c(median = median(x), 
                      l95 = quantile(x, 0.025), 
                      u95 = quantile(x, 0.975))

obs_summary <- 100 * summ(obs_mean_all)  # %
cf_summary  <- 100 * summ(cf_mean_all)   # %
diff_summary <- 100 * summ(100 * diff_dur)  # percentage points

list(obs = obs_summary, cf = cf_summary, diff = diff_summary)


# Keep only people who had DurMin > 0 in the factual world
data2 <- data %>%
  filter(DurMin > 0)

# Observed predictions for that subset
obs_pred_dur2 <- posterior_epred(
  model_observed_dur,
  newdata = data2,
  re_formula = NA
)

# Counterfactual: same people but with duration set to 0
cf_data2 <- data2 %>%
  mutate(DurMin = 0)

cf_pred_dur2 <- posterior_epred(
  model_observed_dur,
  newdata = cf_data2,
  re_formula = NA
)


obs_by_freq2 <- group_means(obs_pred_dur2, data2$travel_frequency)
cf_by_freq2  <- group_means(cf_pred_dur2,  data2$travel_frequency)
diff_by_freq2 <- cf_by_freq2 - obs_by_freq2
cate_table_contact_only <- t(apply(diff_by_freq2, 2, function(col) summ(col)))
cate_table_contact_only

# Mean per draw for this subset
obs_mean_all2 <- rowMeans(obs_pred_dur2)
cf_mean_all2  <- rowMeans(cf_pred_dur2)

# Difference in probability (CF − Obs)
diff_dur2 <- cf_mean_all2 - obs_mean_all2
quantile(diff_dur2, probs = c(0.025, 0.5, 0.975))

obs_summary <- 100 * summ(obs_mean_all2)  # %
cf_summary  <- 100 * summ(cf_mean_all2)   # %
diff_summary <- 100 * summ(100 * diff_dur2)  # percentage points

list(obs = obs_summary, cf = cf_summary, diff = diff_summary)


library(ggplot2)
library(tibble)

# --- Histogram for population-wide effect ---
df_all <- tibble(
  prob = c(obs_mean_all, cf_mean_all),
  type = rep(c("Observed", "Counterfactual"), each = length(obs_mean_all))
)

p_all3 <-ggplot(df_all, aes(x = prob, fill = type)) +
  geom_histogram(position = "identity", alpha = 0.6, bins = 30, colour = "black") +
  scale_fill_manual(values = c("Observed" = "#3B4BA3", "Counterfactual" = "#F4A261")) +
  labs(
    title = "Predicted Infection Probability (Whole Population)",
    x = "Mean infection probability",
    y = "Number of posterior draws"
  ) +
  theme_minimal(base_size = 20)+ theme(
    text = element_text(size = 20, face = "bold"),
    plot.background = element_rect(fill = "white", colour = NA)
  )

# --- Histogram for subset with DurMin > 0 ---
df_subset <- tibble(
  prob = c(obs_mean_all2, cf_mean_all2),
  type = rep(c("Observed", "Counterfactual"), each = length(obs_mean_all2))
)

p_subset <- ggplot(df_subset, aes(x = prob, fill = type)) +
  geom_histogram(position = "identity", alpha = 0.6, bins = 30, colour = "black") +
  scale_fill_manual(values = c("Observed" = "#3B4BA3", "Counterfactual" = "#F4A261")) +
  labs(
    title = "Predicted Infection Probability (Only Those with DurMin > 0)",
    x = "Mean infection probability",
    y = "Number of posterior draws"
  ) +
  theme_minimal(base_size = 14) +
  theme_minimal(base_size = 20)+ theme(
    text = element_text(size = 20, face = "bold"),
    plot.background = element_rect(fill = "white", colour = NA)
  )

# Print the plots
p_all3
p_subset


ggsave("population_effect_duration.png", p_all3, width = 8, height = 6, dpi = 300)
ggsave("subset_effect_duration.png", p_subset, width = 8, height = 6, dpi = 300)

