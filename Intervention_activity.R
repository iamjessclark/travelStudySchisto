rm(list = ls())

# --- Libraries
library(tidyverse)
library(brms)

# --- Data
data <- read.csv("Travel_prevEDIT.csv")
nrow(data)  # 586

# Harmonise types used by the model
data$travel_frequency <- as.factor(data$travel_frequency)
if (!is.factor(data$ActNEW)) data$ActNEW <- factor(data$ActNEW)

# --- Priors
priors <- get_custom_prior("bernoulli")

# --- Fit observed model (Activity as exposure)
model_observed_act <- brm(
  infected_bin ~ ActNEW + age_class + Location + Sex_bin + travel_frequency,
  data = data,
  family = bernoulli(),
  prior = priors,
  iter = 10000, warmup = 3000, chains = 4, cores = 4,
  control = list(max_treedepth = 15)
)

# --- Flag the subgroup that had Occupation in the factual world
data <- data %>% mutate(was_occ = ActNEW == "Occupation")

# --- Build counterfactual dataset: Occupation -> None and DurMin -> 0 for those changed
act_levels <- levels(data$ActNEW)
cf_data_act <- data %>%
  mutate(
    ActNEW = if_else(was_occ, "None", as.character(ActNEW)),
    DurMin = if_else(was_occ, 0, DurMin)  # downstream mediator set to 0 for those changed
  ) %>%
  mutate(
    ActNEW = factor(ActNEW, levels = act_levels)  # preserve original factor levels
  )

# --- Posterior predictions (draws x individuals)
obs_all_draws <- posterior_epred(model_observed_act, newdata = data,       re_formula = NA)
cf_all_draws  <- posterior_epred(model_observed_act, newdata = cf_data_act, re_formula = NA)

# --- Summaries
summ <- function(x){
  x <- as.numeric(x); x <- x[is.finite(x)]
  c(median = stats::median(x),
    l95    = stats::quantile(x, 0.025, names = FALSE),
    u95    = stats::quantile(x, 0.975, names = FALSE))
}

# Whole population (g-computation: mean over individuals per draw)
obs_mean_all <- rowMeans(obs_all_draws)
cf_mean_all  <- rowMeans(cf_all_draws)

pop_obs_summary  <- 100 * summ(obs_mean_all)
pop_cf_summary   <- 100 * summ(cf_mean_all)
pop_diff_summary <- 100 * summ(cf_mean_all - obs_mean_all)  # percentage points

# Subgroup: those who were Occupation in the factual world
was_occ <- data$was_occ
stopifnot(any(was_occ))

obs_occ_draws <- rowMeans(obs_all_draws[, was_occ, drop = FALSE], na.rm = TRUE)
cf_occ_draws  <- rowMeans(cf_all_draws[,  was_occ, drop = FALSE], na.rm = TRUE)

occ_obs_summary  <- 100 * summ(obs_occ_draws)
occ_cf_summary   <- 100 * summ(cf_occ_draws)
occ_diff_summary <- 100 * summ(cf_occ_draws - obs_occ_draws)

# --- Table for manuscript text
results_tbl <- tibble::tibble(
  group      = rep(c("Whole population", "Occupation subgroup"), each = 3),
  scenario   = rep(c("Observed", "Counterfactual", "Difference (CF−Obs)"), times = 2),
  median_pct = c(pop_obs_summary["median"],  pop_cf_summary["median"],  pop_diff_summary["median"],
                 occ_obs_summary["median"],  occ_cf_summary["median"],  occ_diff_summary["median"]),
  l95_pct    = c(pop_obs_summary["l95"],     pop_cf_summary["l95"],     pop_diff_summary["l95"],
                 occ_obs_summary["l95"],     occ_cf_summary["l95"],     occ_diff_summary["l95"]),
  u95_pct    = c(pop_obs_summary["u95"],     pop_cf_summary["u95"],     pop_diff_summary["u95"],
                 occ_obs_summary["u95"],     occ_cf_summary["u95"],     occ_diff_summary["u95"])
)
print(results_tbl)

# --- Plots (same style as your travel intervention)
library(ggplot2)

# Whole population
df_all <- tibble(
  prob = c(obs_mean_all, cf_mean_all),
  type = rep(c("Observed", "Counterfactual"), each = length(obs_mean_all))
)

p_all <- ggplot(df_all, aes(x = prob, fill = type)) +
  geom_histogram(position = "identity", alpha = 0.6, bins = 30, colour = "black") +
  scale_fill_manual(values = c("Observed" = "#3B4BA3", "Counterfactual" = "#F4A261")) +
  labs(
    title = "Predicted Infection Probability (Whole Population)\nIntervention: Occupation → None",
    x = "Mean infection probability",
    y = "Number of posterior draws"
  ) +
  theme_minimal(base_size = 20)+ theme(
    text = element_text(size = 20, face = "bold"),
    plot.background = element_rect(fill = "white", colour = NA)
  )

# Occupation subgroup
df_occ <- tibble(
  prob = c(obs_occ_draws, cf_occ_draws),
  type = rep(c("Observed", "Counterfactual"), each = length(obs_occ_draws))
)

p_occ <- ggplot(df_occ, aes(x = prob, fill = type)) +
  geom_histogram(position = "identity", alpha = 0.6, bins = 30, colour = "black") +
  scale_fill_manual(values = c("Observed" = "#3B4BA3", "Counterfactual" = "#F4A261")) +
  labs(
    title = "Predicted Infection Probability (Occupation Subgroup)\nIntervention: Occupation → None",
    x = "Mean infection probability",
    y = "Number of posterior draws"
  ) +
  theme_minimal(base_size = 20)+ theme(
    text = element_text(size = 20, face = "bold"),
    plot.background = element_rect(fill = "white", colour = NA)
  )

# Draw plots
p_all
p_occ


ggsave("population_effect_act.png", p_all, width = 8, height = 6, dpi = 300)
ggsave("subset_effect_act.png", p_occ, width = 8, height = 6, dpi = 300)

###worm buden and domestic
###
filt <- data%>%
  filter(mean_sm >0)

Act_mean_sm <- brm(mean_sm ~ ActNEW+ age_class+Location+Sex_bin+travel_frequency,
                   data=filt,
                   family="gamma",
                   prior=  c(
                     prior(normal(2.5, 1), class = "Intercept"),
                     prior(normal(0, 10), class = "b"),
                     prior(exponential(1), class = "shape")),
                   iter = 10000, warmup = 3000, chains = 4, cores = 8,
                   control = list(max_treedepth = 15))


# --- Flag Domestic cases in the factual data
filt <- filt %>%
  mutate(was_domestic = ActNEW == "Domestic")

# Store factor levels
act_levels <- levels(filt$ActNEW)

# --- Build counterfactual data (Domestic → None)
cf_data_dom <- filt %>%
  mutate(
    ActNEW = if_else(was_domestic, "None", as.character(ActNEW)),
    DurMin = if_else(was_domestic, 0, DurMin)  # downstream duration set to 0
  ) %>%
  mutate(
    ActNEW = factor(ActNEW, levels = act_levels)
  )

# --- Posterior predictions (draws x individuals)
obs_all_draws <- posterior_epred(Act_mean_sm, newdata = filt,       re_formula = NA)
cf_all_draws  <- posterior_epred(Act_mean_sm, newdata = cf_data_dom, re_formula = NA)


# 2. Mask extreme values (>150 EPG) with NA
obs_all_draws[obs_all_draws > 150] <- NA
cf_all_draws[cf_all_draws > 150]   <- NA


# --- Summaries
summ <- function(x){
  x <- as.numeric(x)
  x <- x[is.finite(x)]                # drop NA/NaN/Inf
  stopifnot(length(x) > 0)            # fail fast if empty
  c(
    median = stats::median(x),
    l95    = stats::quantile(x, 0.025, names = FALSE),
    u95    = stats::quantile(x, 0.975, names = FALSE)
  )
}

# Whole population
obs_mean_all <- rowMeans(obs_all_draws)
cf_mean_all  <- rowMeans(cf_all_draws)

pop_obs_summary  <- summ(obs_mean_all)
pop_cf_summary   <- summ(cf_mean_all)
pop_diff_summary <- summ(cf_mean_all - obs_mean_all)  # absolute difference in mean EPG

# Domestic subgroup
was_domestic <- filt$was_domestic
obs_dom_draws <- rowMeans(obs_all_draws[, was_domestic, drop = FALSE], na.rm = TRUE)
cf_dom_draws  <- rowMeans(cf_all_draws[,  was_domestic, drop = FALSE], na.rm = TRUE)

dom_obs_summary  <- summ(obs_dom_draws)
dom_cf_summary   <- summ(cf_dom_draws)
dom_diff_summary <- summ(cf_dom_draws - obs_dom_draws)

# --- Results table
results_tbl <- tibble::tibble(
  group      = rep(c("Whole population", "Domestic subgroup"), each = 3),
  scenario   = rep(c("Observed", "Counterfactual", "Difference (CF−Obs)"), times = 2),
  median     = c(pop_obs_summary["median"],  pop_cf_summary["median"],  pop_diff_summary["median"],
                 dom_obs_summary["median"],  dom_cf_summary["median"],  dom_diff_summary["median"]),
  l95        = c(pop_obs_summary["l95"],     pop_cf_summary["l95"],     pop_diff_summary["l95"],
                 dom_obs_summary["l95"],     dom_cf_summary["l95"],     dom_diff_summary["l95"]),
  u95        = c(pop_obs_summary["u95"],     pop_cf_summary["u95"],     pop_diff_summary["u95"],
                 dom_obs_summary["u95"],     dom_cf_summary["u95"],     dom_diff_summary["u95"])
)
print(results_tbl)

# --- Plots
library(ggplot2)

# Whole population
df_all <- tibble(
  value = c(obs_mean_all, cf_mean_all),
  type  = rep(c("Observed", "Counterfactual"), each = length(obs_mean_all))
)

p_all2 <- ggplot(df_all, aes(x = value, fill = type)) +
  geom_histogram(position = "identity", alpha = 0.6, bins = 30, colour = "black") +
  scale_fill_manual(values = c("Observed" = "#3B4BA3", "Counterfactual" = "#F4A261")) +
  labs(
    title = "Predicted Mean EPG (Whole Population)\nIntervention: Domestic → None",
    x = "Mean EPG",
    y = "Number of posterior draws"
  ) +
  theme_minimal(base_size = 20)+ theme(
    text = element_text(size = 20, face = "bold"),
    plot.background = element_rect(fill = "white", colour = NA)
  )

# Domestic subgroup
df_dom <- tibble(
  value = c(obs_dom_draws, cf_dom_draws),
  type  = rep(c("Observed", "Counterfactual"), each = length(obs_dom_draws))
)

p_dom <- ggplot(df_dom, aes(x = value, fill = type)) +
  geom_histogram(position = "identity", alpha = 0.6, bins = 30, colour = "black") +
  scale_fill_manual(values = c("Observed" = "#3B4BA3", "Counterfactual" = "#F4A261")) +
  labs(
    title = "Predicted Mean EPG (Domestic Subgroup)\nIntervention: Domestic → None",
    x = "Mean EPG",
    y = "Number of posterior draws"
  ) +
  theme_minimal(base_size = 20)+ theme(
    text = element_text(size = 20, face = "bold"),
    plot.background = element_rect(fill = "white", colour = NA)
  )

p_all2
p_dom
ggsave("population_effect_EPG_act.png", p_all2, width = 8, height = 6, dpi = 300)
ggsave("subset_effect_EPG_act.png", p_dom, width = 8, height = 6, dpi = 300)
