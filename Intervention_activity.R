# Load required libraries
library(tidyverse)
library(brms)

#=================================================================
#ACTIVITY INTERVENTION ########################
#================================================================
# Set priors
priors <- get_custom_prior("bernoulli")

#

# 1. Fit observed model (total effect with covariates)
model_observed_act <- brm(
  infected_bin ~ActNEW+age_class+Location+Sex_bin+travel_frequency,
  data = data,
  family = bernoulli(),
  prior = priors,
  iter = 10000, warmup = 3000, chains = 4, cores = 4,
  control = list(max_treedepth = 15)
)

# 2. Fit counterfactual model (total effect, backdoor paths are now closed so no confounders)
model_cf_act<- brm(
  infected_bin ~ ActNEW,
  data = data,
  family = bernoulli(),
  prior = priors,
  iter = 10000, warmup = 3000, chains = 4, cores = 4,
  control = list(max_treedepth = 15)
)


# # 3. Predict infection under observed model (original frequency)
# obs_pred_act <- rowMeans(posterior_epred(
#   model_observed_act,
#   newdata = data,
#   re_formula = NA
# ))

#new step 3
## Posterior predictions: full uncertainty (individual outcome variability)
obs_ppred_act <- posterior_predict(
  model_observed_act,
  newdata = data,
  re_formula = NA
)

# Average predicted infection *per individual* across posterior draws
obs_ppred_act_mean <- rowMeans(obs_ppred_act)

# 4. Create counterfactual: Everyone doing occupation is stopped
data_cf_act_OCC <- data %>%
  mutate(ActNEW = case_when(
    ActNEW == "Occupation" ~ "None",
    TRUE ~ ActNEW
  ))

# 
# # 5. Predict infection under counterfactual model
# cf_pred_act <- rowMeans(posterior_epred(
#   model_cf_act,
#   newdata = data_cf_act_OCC,
#   re_formula = NA
# ))

#new step 5
cf_ppred_act <- posterior_predict(
  model_cf_act,
  newdata = data_cf_act_OCC,
  re_formula = NA
)

cf_ppred_act_mean <- rowMeans(cf_ppred_act)


# # 6. Summarise change in prevalence 
# diff_act<- cf_pred_act - obs_pred_act
# quantile(diff_dur, probs = c(0.025, 0.5, 0.975))


#new step 6
diff_ppred_act <- cf_ppred_act_mean - obs_ppred_act_mean
quantile(diff_ppred_act, probs = c(0.025, 0.5, 0.975))


########################### activty ==none


# 4. Create counterfactual: no duration in the lake for everyone
data_cf_act_OC2C <- data %>%
  mutate(ActNEW = case_when(
    ActNEW !="None" ~ "None",
    TRUE ~ ActNEW
  ))

# 5. Predict infection under counterfactual model
cf_pred_act2 <- rowMeans(posterior_epred(
  model_cf_act,
  newdata = data_cf_act_OC2C,
  re_formula = NA
))

# 6. Summarise change in prevalence 
diff_act2<- cf_pred_act2 - obs_pred_act
quantile(diff_act2, probs = c(0.025, 0.5, 0.975))

#check by just comparing coefficents for each activity category

posterior <- as_draws_df(model_cf_act)

occupation_effect <- posterior$b_ActNEWOccupation

occupation_prob_diff <- plogis(posterior$b_Intercept + occupation_effect) - 
  plogis(posterior$b_Intercept)

quantile(occupation_prob_diff, c(0.025, 0.5, 0.975))

#====================================================================================
# Graphs     #######
# ===================================================================================
#
library(tidyverse)

####
#Occupation to none

# Prepare data
df_act_Occ <- tibble(
  prob = c(obs_pred_act, cf_pred_act),
  type = rep(c("Observed", "Counterfactual"), each = length(obs_pred_act))
)

# Plot
ggplot(df_act_Occ, aes(x = prob, fill = type)) +
  geom_histogram(position = "identity", alpha = 0.6, bins = 30, colour = "black") +
  scale_fill_manual(values = c("Observed" = "#3B4BA3", "Counterfactual" = "#F4A261")) +
  labs(
    title = "Predicted Infection Probability (Occupation changed to none)",
    x = "Infection probability",
    y = "Number of posterior predictions"
  ) +
  theme_minimal()+
  theme(text=element_text(size=18),
        axis.title=element_text(face="bold"))

#Histogram for everyone reporting no activity

# Prepare data
df_act <- tibble(
  prob = c(obs_pred_act, cf_pred_act2),
  type = rep(c("Observed", "Counterfactual"), each = length(obs_pred_act))
)

# Plot
ggplot(df_act, aes(x = prob, fill = type)) +
  geom_histogram(position = "identity", alpha = 0.6, bins = 30, colour = "black") +
  scale_fill_manual(values = c("Observed" = "#3B4BA3", "Counterfactual" = "#F4A261")) +
  labs(
    title = "Predicted Infection Probability (No activity)",
    x = "Infection probability",
    y = "Number of posterior predictions"
  ) +
  theme_minimal()+
  theme(text=element_text(size=18),
        axis.title=element_text(face="bold"))

#side by side graphs
#
#
#library(tidyverse)
library(patchwork)  # for side-by-side plots

# Simulated draws (from your posterior_epred output)
# Make sure you have these ready:
# - obs_pred_all, cf_pred_all: full population (n_draws × n_individuals)
# - obs_pred_daily, cf_pred_daily: daily travellers (n_draws × n_individuals)
# 
# 
# 

# Convert to long format
make_long_df_act <- function(pred_matrix, label) {
  as.data.frame(pred_matrix) %>%
    pivot_longer(everything(), names_to = "Individual", values_to = "Infection_Prob") %>%
    mutate(Type = label)
}

df_all_obs_act <- make_long_df(obs_pred_act, "Observed")
df_all_cf_act  <- make_long_df(cf_pred_act,  "Counterfactual")

df_act2_cf  <- make_long_df(cf_pred_act2,  "Counterfactual")

# Combine
df_all_combined_actOcc  <- bind_rows(df_all_obs_act, df_all_cf_act)
df_act2_combinedALL <- bind_rows(df_all_obs_act, df_act2_cf)

# Plot settings
plot_hist2 <- function(df, title) {
  ggplot(df, aes(x = Infection_Prob, fill = Type)) +
    geom_histogram(position = "identity", alpha = 0.6, bins = 30, colour = "black") +
    scale_fill_manual(name = NULL,
                      values = c("Observed" = "#4477AA", "Counterfactual" = "#EE7733")) +
    labs(
      x = "",
      y = "",
      title = title,
      
    ) +
    theme_minimal(base_size = 13)+
    theme(text=element_text(size=28, face="bold"), 
          axis.title = element_text(face="bold"))
}

p2 <- plot_hist2(df_all_combined_actOcc ,   "Occupation only")
p1 <- plot_hist2(df_act2_combinedALL, "All activities")

# Side-by-side layout
side_act <-p1 + p2 + plot_layout(guides = "collect") & theme(legend.position = "bottom")

side_act

ggsave("Figures/Activity_intervention.png", side_act, width = 12, height = 10, dpi = 300)

#################
################## 
#What about on worm burden

