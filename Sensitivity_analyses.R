# Sensitivity_analyses.R
# Robustness analyses addressing four peer reviewer comments
#
# Comment 1: Cross-sectional design -- no temporal ordering, unmeasured confounders
#            -> E-value analysis (unmeasured confounding threshold)
#
# Comment 2: Self-reported exposures prone to social-desirability bias; MDA
#            participation may introduce residual confounding
#            -> Measurement error bounding for DurMin (classical attenuation formula)
#            -> E-value treating MDA participation as potential confounder
#            -> Travel effect with vs without MDA in adjustment set
#
# Comment 3: Counterfactual simulations assume independent manipulation; no
#            spillover or behavioural adaptation
#            -> Spillover bounding (range of assumed spillover fractions)
#            -> Behavioural substitution bounding
#
# Comment 4: Activity may not be a true mediator; conditioning on it could remove
#            real travel effects; DAG not empirically verifiable
#            -> Compare travel effect under two DAG structures:
#               Model A = Act is mediator (total effect, do NOT adjust for Act)
#               Model B = Act is confounder/common cause (adjust for Act)
#            -> Present unadjusted travel association alongside both

rm(list = ls())
set.seed(2024)

library(brms)
library(dplyr)
library(ggplot2)
library(tidybayes)
library(tibble)
library(tidyr)
library(posterior)
library(scales)

# Install EValue if needed: install.packages("EValue")
if (!requireNamespace("EValue", quietly = TRUE)) {
  message("Install EValue package: install.packages('EValue')")
  HAS_EVALUE <- FALSE
} else {
  library(EValue)
  HAS_EVALUE <- TRUE
}

# ---- Data preparation -------------------------------------------------------
data <- read.csv("Travel_prevEDIT.csv")
nrow(data)  # 586 (one duplicate)

data$travel_frequency <- as.factor(data$travel_frequency)
data$ActNEW           <- as.factor(data$ActNEW)
data$ActNEW           <- relevel(data$ActNEW, ref = "None")
data$Location         <- as.factor(data$Location)
data$Location         <- relevel(data$Location, ref = "Naawa")
data$infected_bin     <- as.numeric(as.character(data$infected_bin))

# Check levels so parameter names can be verified
cat("travel_frequency levels:", levels(data$travel_frequency), "\n")
cat("ActNEW levels:          ", levels(data$ActNEW), "\n")
cat("Location levels:        ", levels(data$Location), "\n")

# Priors matching priors_weak.R used in original manuscript models
# Using the same get_custom_prior() priors ensures E-value ORs/MRs match
# the published results exactly.
prior_bern <- c(
  prior(normal(0, 10),   class = "b"),
  prior(normal(-1.5, 2.5), class = "Intercept")
)

prior_bern_dur <- c(
  prior(normal(0, 10),   class = "b"),
  prior(normal(-1.5, 2.5), class = "Intercept")
)

prior_gamma <- c(
  prior(normal(0, 10),   class = "b"),
  prior(normal(2.5, 1),  class = "Intercept"),
  prior(exponential(1),  class = "shape")
)

# ---- Helper: save / load brms models ----------------------------------------
model_dir <- "models/sensitivity"
if (!dir.exists(model_dir)) dir.create(model_dir, recursive = TRUE)

load_or_fit <- function(filename, formula, data, family, prior, iter = 10000,
                         warmup = 3000, chains = 4, cores = 4,
                         control = list(max_treedepth = 15)) {
  path <- file.path(model_dir, filename)
  if (file.exists(path)) {
    message("Loading cached model: ", filename)
    return(readRDS(path))
  }
  message("Fitting model: ", filename, " (this may take several minutes)")
  m <- brm(formula = formula, data = data, family = family, prior = prior,
            iter = iter, warmup = warmup, chains = chains, cores = cores,
            control = control, seed = 2024)
  saveRDS(m, path)
  m
}

# ---- Helper: extract posterior OR for a named coefficient -------------------
posterior_or <- function(model, coef_name) {
  d <- as_draws_df(model)
  if (!coef_name %in% names(d)) {
    stop("Coefficient '", coef_name, "' not found. Available: ",
         paste(grep("^b_", names(d), value = TRUE), collapse = ", "))
  }
  vals <- d[[coef_name]]
  tibble(
    log_or_med = median(vals),
    log_or_lo  = as.numeric(quantile(vals, 0.025)),
    log_or_hi  = as.numeric(quantile(vals, 0.975)),
    or_med     = exp(median(vals)),
    or_lo      = exp(as.numeric(quantile(vals, 0.025))),
    or_hi      = exp(as.numeric(quantile(vals, 0.975)))
  )
}

# ---- Manual E-value formula -------------------------------------------------
# E = RR + sqrt(RR * (RR - 1))
# For OR with non-rare outcome: approximate RR from OR and baseline risk P0
evalue_from_or <- function(or, or_lo, p0) {
  # Convert OR to approximate RR using baseline risk (P0 = P(outcome | unexposed))
  rr_from_or <- function(or, p0) or / ((1 - p0) + p0 * or)
  ev <- function(rr) {
    if (rr < 1) rr <- 1 / rr  # flip for protective
    rr + sqrt(rr * (rr - 1))
  }
  rr_est  <- rr_from_or(or, p0)
  rr_lb   <- rr_from_or(or_lo, p0)
  ev_est  <- ev(rr_est)
  ev_lb   <- if (or_lo > 1) ev(rr_lb) else 1  # E-value for lower bound
  list(rr_estimate  = rr_est,
       rr_lower_ci  = rr_lb,
       evalue_est   = ev_est,
       evalue_lo_ci = ev_lb,
       p0_used      = p0)
}

# =============================================================================
# SENSITIVITY 1: COMPREHENSIVE E-VALUE TABLE — ALL KEY ASSOCIATIONS
# Covers all main results: travel (daily), activity (all types), MDA
# Outcomes: infection prevalence (binary) and burden (EPG, Gamma, infected only)
# Adjustment sets from DAG:
#   Travel -> Inf/Burd:  {Age, Sex, Location}
#   Act    -> Inf/Burd:  {Age, Sex, Location}
#   MDA    -> Inf/Burd:  {Travel_freq, Age, Location}  [Sex not in adj. set for MDA]
# =============================================================================
message("\n=== SENSITIVITY 1: COMPREHENSIVE E-VALUE TABLE ===\n")

# ---- Filtered data for burden models (infected individuals only) ----
filt <- data %>% filter(mean_sm > 0)

# ---- Fit all 6 models -------------------------------------------------------
model_travel_inf <- load_or_fit(
  "model_travel_inf.rds",
  formula = infected_bin ~ travel_frequency + age_class + Sex_bin + Location,
  data = data, family = bernoulli(), prior = prior_bern
)

model_act_inf <- load_or_fit(
  "model_act_inf.rds",
  formula = infected_bin ~ ActNEW + age_class + Sex_bin + Location + travel_frequency,
  data = data, family = bernoulli(), prior = prior_bern
)

model_mda_inf <- load_or_fit(
  "model_mda_inf.rds",
  formula = infected_bin ~ MDA + travel_frequency + age_class + Location,
  data = data, family = bernoulli(), prior = prior_bern
)

model_travel_burd <- load_or_fit(
  "model_travel_burd.rds",
  formula = mean_sm ~ travel_frequency + age_class + Sex_bin + Location,
  data = filt, family = Gamma(link = "log"), prior = prior_gamma
)

model_act_burd <- load_or_fit(
  "model_act_burd.rds",
  formula = mean_sm ~ ActNEW + age_class + Sex_bin + Location + travel_frequency,
  data = filt, family = Gamma(link = "log"), prior = prior_gamma
)

model_mda_burd <- load_or_fit(
  "model_mda_burd.rds",
  formula = mean_sm ~ MDA + travel_frequency + age_class + Location,
  data = filt, family = Gamma(link = "log"), prior = prior_gamma
)

# Keep alias so Sensitivities 2-4 (which use model_travel_tot) still work
model_travel_tot <- model_travel_inf
model_act_tot    <- model_act_inf

# ---- Baseline risks for OR-to-RR conversion (binary models only) ------------
p0_never  <- mean(data$infected_bin[data$travel_frequency == "0"], na.rm = TRUE)
p0_no_act <- mean(data$infected_bin[data$ActNEW == "None"],        na.rm = TRUE)
p0_no_mda <- mean(data$infected_bin[data$MDA == 0],                na.rm = TRUE)

cat("Baseline risks — never travel:", round(p0_never * 100, 1),
    "%, no activity:", round(p0_no_act * 100, 1),
    "%, no MDA:", round(p0_no_mda * 100, 1), "%\n\n")

# ---- General E-value function -----------------------------------------------
# type = "binary": exp(coef) is an OR; convert to approx RR via baseline risk p0
# type = "gamma":  exp(coef) is a mean ratio; used directly (same formula as RR)
# The "inner CI bound" is the CI bound closest to null:
#   harmful (effect > 1): lower CI;  protective (effect < 1): upper CI
compute_evalue_row <- function(exposure, outcome_label, model, coef_name,
                                type = c("binary", "gamma"), p0 = NULL) {
  type <- match.arg(type)
  d    <- as_draws_df(model)
  if (!coef_name %in% names(d)) {
    warning("Coefficient not found: ", coef_name, " in model for ", outcome_label)
    return(NULL)
  }
  vals    <- d[[coef_name]]
  log_med <- median(vals)
  log_lo  <- as.numeric(quantile(vals, 0.025))
  log_hi  <- as.numeric(quantile(vals, 0.975))

  eff_med <- exp(log_med); eff_lo <- exp(log_lo); eff_hi <- exp(log_hi)

  if (type == "binary") {
    to_rr  <- function(or) or / ((1 - p0) + p0 * or)
    rr_med <- to_rr(eff_med); rr_lo <- to_rr(eff_lo); rr_hi <- to_rr(eff_hi)
    eff_label <- "OR"
  } else {
    rr_med <- eff_med; rr_lo <- eff_lo; rr_hi <- eff_hi
    eff_label <- "MR"
  }

  ev_fn <- function(rr) {
    if (!is.finite(rr) || rr == 1) return(1)
    if (rr < 1) rr <- 1 / rr
    rr + sqrt(rr * (rr - 1))
  }

  # CI bound closest to null
  if (rr_med >= 1) {
    inner_rr <- rr_lo; cri_clears_null <- rr_lo > 1
  } else {
    inner_rr <- rr_hi; cri_clears_null <- rr_hi < 1
  }

  tibble(
    exposure          = exposure,
    outcome           = outcome_label,
    effect_type       = eff_label,
    effect_est        = round(eff_med, 2),
    cri_lo            = round(eff_lo,  2),
    cri_hi            = round(eff_hi,  2),
    rr_or_mr_approx   = round(rr_med,  2),
    evalue_estimate   = round(ev_fn(rr_med), 2),
    evalue_inner_cri  = round(if (cri_clears_null) ev_fn(inner_rr) else 1, 2)
  )
}

# Applies compute_evalue_row to every coefficient matching a variable prefix
extract_evalue_rows <- function(var_prefix, label_fn, model, type, p0 = NULL,
                                 outcome_label) {
  d     <- as_draws_df(model)
  coefs <- grep(paste0("^b_", var_prefix), names(d), value = TRUE)
  bind_rows(lapply(coefs, function(coef) {
    level <- sub(paste0("b_", var_prefix), "", coef)
    compute_evalue_row(label_fn(level), outcome_label, model, coef, type, p0)
  }))
}

# ---- Build the comprehensive supplementary table ----------------------------
daily_coef <- "b_travel_frequency92"

evalue_supp <- bind_rows(

  # ---- INFECTION (prevalence, binary) ----------------------------------------
  compute_evalue_row(
    "Daily travel (92×/3mo) vs. never", "Infection (prevalence)",
    model_travel_inf, daily_coef, "binary", p0_never),

  extract_evalue_rows(
    "ActNEW", function(l) paste0(l, " vs. no lake activity"),
    model_act_inf, "binary", p0_no_act, "Infection (prevalence)"),

  extract_evalue_rows(
    "MDA", function(l) "MDA participation vs. none",
    model_mda_inf, "binary", p0_no_mda, "Infection (prevalence)"),

  # ---- BURDEN (mean EPG, Gamma, infected only) --------------------------------
  compute_evalue_row(
    "Daily travel (92×/3mo) vs. never", "Burden (EPG, infected only)",
    model_travel_burd, daily_coef, "gamma"),

  extract_evalue_rows(
    "ActNEW", function(l) paste0(l, " vs. no lake activity"),
    model_act_burd, "gamma", NULL, "Burden (EPG, infected only)"),

  extract_evalue_rows(
    "MDA", function(l) "MDA participation vs. none",
    model_mda_burd, "gamma", NULL, "Burden (EPG, infected only)")
)

cat("=== Supplementary E-value table ===\n")
cat("effect_type: OR = odds ratio (binary models); MR = mean ratio (Gamma models)\n")
cat("evalue_inner_cri: E-value for CI bound closest to null\n",
    "  (lower CI for harmful effects, upper CI for protective effects)\n\n")
print(evalue_supp, n = Inf)

write.csv(evalue_supp, "Figures/evalue_supplementary_table.csv", row.names = FALSE)
cat("Supplementary table saved: Figures/evalue_supplementary_table.csv\n\n")

# EValue package cross-check on the travel result if available
if (HAS_EVALUE) {
  travel_row <- evalue_supp %>%
    filter(exposure == "Daily travel (92×/3mo) vs. never",
           outcome  == "Infection (prevalence)")
  cat("EValue package cross-check (travel -> infection):\n")
  print(evalues.OR(est  = travel_row$effect_est,
                   lo   = travel_row$cri_lo,
                   hi   = travel_row$cri_hi,
                   rare = FALSE))
}


# =============================================================================
# SENSITIVITY 2a: MEASUREMENT ERROR IN SELF-REPORTED WATER CONTACT DURATION
# Classical measurement error in a continuous predictor attenuates the
# coefficient toward zero. True effect = observed / reliability_ratio.
# =============================================================================
message("\n=== SENSITIVITY 2a: MEASUREMENT ERROR IN DurMin ===\n")

# Adjustment set for Dur -> Inf (total & direct same): Act, Age, Sex, Location
model_dur_tot <- load_or_fit(
  "model_dur_total.rds",
  formula  = infected_bin ~ DurMin + ActNEW + age_class + Sex_bin + Location,
  data     = data,
  family   = bernoulli(),
  prior    = prior_bern_dur
)

summary(model_dur_tot)

draws_dur <- as_draws_df(model_dur_tot)$b_DurMin
beta_dur  <- c(med = median(draws_dur), lo = quantile(draws_dur, 0.025),
               hi  = quantile(draws_dur, 0.975))
cat("Observed DurMin log-OR per minute: median =", round(beta_dur["med"], 4),
    " 95% CrI:", round(beta_dur["lo"], 4), "to", round(beta_dur["hi"], 4), "\n\n")

# Classical measurement error: observed beta = lambda * true beta
# lambda = reliability ratio = Var_true / (Var_true + Var_error)
# If self-reports are unreliable, lambda < 1 and true effect is LARGER.
# This bias is toward the null -> our conclusions are conservative.

lambda_vals <- c(0.5, 0.6, 0.7, 0.8, 0.9, 1.0)

me_sensitivity <- tibble(
  reliability_ratio     = lambda_vals,
  pct_measurement_error = paste0(round((1 - lambda_vals) * 100), "%"),
  log_or_corrected_med  = beta_dur["med"] / lambda_vals,
  log_or_corrected_lo   = beta_dur["lo"]  / lambda_vals,
  log_or_corrected_hi   = beta_dur["hi"]  / lambda_vals,
  or_per_min_corrected  = exp(beta_dur["med"] / lambda_vals),
  direction             = "Attenuation: true effect at least as large as corrected value"
)

cat("=== Measurement error sensitivity for DurMin ===\n")
cat("(lambda = reliability ratio; 1.0 = perfect measurement)\n\n")
print(me_sensitivity %>% select(-direction), n = Inf)
write.csv(me_sensitivity, "Figures/measurement_error_sensitivity.csv", row.names = FALSE)


# Sensitivity 2b has been merged into Sensitivity 4 (DAG structural uncertainty).
# All mediator-as-confounder tests (MDA, Activity, Duration, and combinations)
# are now handled together in a single comprehensive structural sensitivity analysis.


# =============================================================================
# SENSITIVITY 3: SPILLOVER / BEHAVIOURAL ADAPTATION BOUNDING
# Counterfactual: all daily travellers (92×) set to zero travel
# Base case assumes: (a) no spillover to non-travellers, (b) no substitution
# Sensitivity: bound estimates under range of spillover / substitution fracs
# =============================================================================
message("\n=== SENSITIVITY 3: SPILLOVER BOUNDING ===\n")

# Construct counterfactual data (same as Intervention_travel.R)
data <- data %>% mutate(was_daily = travel_frequency == "92")

tf_levels  <- levels(data$travel_frequency)
act_levels <- if (is.factor(data$ActNEW)) levels(data$ActNEW) else NULL

cf_data_daily <- data %>%
  mutate(
    travel_frequency = if_else(was_daily, "0", as.character(travel_frequency)),
    DurMin  = if_else(travel_frequency == "0", 0, DurMin),
    ActNEW  = if_else(travel_frequency == "0", "None", as.character(ActNEW))
  ) %>%
  mutate(
    travel_frequency = factor(travel_frequency, levels = tf_levels),
    ActNEW = if (!is.null(act_levels)) factor(ActNEW, levels = act_levels) else ActNEW
  )

obs_all  <- posterior_epred(model_travel_tot, newdata = data,          re_formula = NA)
cf_all   <- posterior_epred(model_travel_tot, newdata = cf_data_daily, re_formula = NA)

obs_pop  <- rowMeans(obs_all)
cf_pop   <- rowMeans(cf_all)
diff_pop <- cf_pop - obs_pop  # negative = reduction in prevalence

summ_fn <- function(x) {
  x <- as.numeric(x[is.finite(x)])
  c(med = 100 * median(x), lo = 100 * quantile(x, 0.025, names = FALSE),
    hi  = 100 * quantile(x, 0.975, names = FALSE))
}

base_diff <- summ_fn(diff_pop)
cat("Base counterfactual: if all daily travellers stopped travelling\n")
cat("  Population prevalence difference (CF - obs):",
    round(base_diff["med"], 2), "% [",
    round(base_diff["lo"], 2), ",", round(base_diff["hi"], 2), "]\n\n")

# --- 3a: SPILLOVER to non-daily travellers ---------------------------------
# Scenario: stopping travel in daily travellers also benefits non-travellers
# (e.g., reduced environmental shedding/infectious pressure in shared waterways)
# This makes the intervention MORE effective than the base case.
# Spillover fraction = fraction of the daily-traveller risk reduction that
# also applies to non-daily-travellers.

not_daily <- !data$was_daily
spillover_fracs <- c(0, 0.10, 0.20, 0.30)

# Per-draw mean risk reduction for daily travellers (negative = lower risk in CF)
per_person_reduction    <- cf_all[, data$was_daily, drop = FALSE] -
  obs_all[, data$was_daily, drop = FALSE]
mean_reduction_per_draw <- rowMeans(per_person_reduction)

spillover_results <- bind_rows(lapply(spillover_fracs, function(sf) {
  adj_cf <- cf_all
  adj_cf[, not_daily] <- obs_all[, not_daily] + mean_reduction_per_draw * sf
  adj_diff <- rowMeans(adj_cf) - obs_pop
  s <- summ_fn(adj_diff)
  tibble(spillover_frac = sf,
         scenario       = paste0(sf * 100, "% spillover benefit to non-daily-travellers"),
         diff_med_pct   = s["med"],
         diff_lo_pct    = s["lo"],
         diff_hi_pct    = s["hi"])
}))

cat("=== Spillover bounding: effect on population prevalence ===\n")
print(spillover_results, n = Inf)

# --- 3b: BEHAVIOURAL SUBSTITUTION ------------------------------------------
# Scenario: some daily travellers substitute travel to other (lower-risk) sites
# so the effective exposure reduction is less than 100%.
# Substitution fraction = fraction of exposure that is merely displaced.
# sub_frac = 0: no substitution (full benefit)
# sub_frac = 0.30: 30% of daily travellers merely substitute elsewhere

sub_fracs <- c(0, 0.10, 0.20, 0.30, 0.50)

substitution_results <- bind_rows(lapply(sub_fracs, function(sf) {
  # Effective difference = (1 - sf) * base_difference per posterior draw
  adj_diff <- (1 - sf) * diff_pop
  s <- summ_fn(adj_diff)
  tibble(sub_frac       = sf,
         scenario       = paste0(sf * 100, "% behavioural substitution"),
         diff_med_pct   = s["med"],
         diff_lo_pct    = s["lo"],
         diff_hi_pct    = s["hi"])
}))

cat("\n=== Behavioural substitution bounding ===\n")
print(substitution_results, n = Inf)

# Combine and plot
spillover_combined <- bind_rows(
  mutate(spillover_results, type = "Spillover to non-travellers"),
  mutate(substitution_results, type = "Behavioural substitution",
         spillover_frac = sub_frac) %>% select(-sub_frac)
)

p_spill <- ggplot(spillover_combined,
                  aes(x = spillover_frac, y = diff_med_pct,
                      ymin = diff_lo_pct, ymax = diff_hi_pct,
                      colour = type, fill = type)) +
  geom_ribbon(alpha = 0.20, colour = NA) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  scale_x_continuous(labels = function(x) paste0(x * 100, "%")) +
  labs(x = "Assumed fraction (spillover or substitution)",
       y = "Population prevalence difference, %\n(counterfactual − observed)",
       colour = NULL, fill = NULL,
       title = "Sensitivity 3: spillover and substitution bounding") +
  theme_minimal(base_size = 13) +
  theme(plot.background  = element_rect(fill = "white", colour = NA),
        legend.position  = "bottom")

ggsave("Figures/spillover_bounding.png", p_spill, width = 8, height = 5, dpi = 300)
cat("Figure saved: Figures/spillover_bounding.png\n")

write.csv(spillover_combined, "Figures/spillover_bounding.csv", row.names = FALSE)


# =============================================================================
# SENSITIVITY 4: DAG STRUCTURAL UNCERTAINTY — ALL MEDIATORS AS POTENTIAL CONFOUNDERS
#
# In the DAG, Travel -> {Activity, MDA, Duration} -> Infection.
# All three are mediators: the correct total-effect model does not condition on any.
# However, each could alternatively be a confounder if driven by unmeasured
# common causes of both travel frequency and infection.
#
# We test all 8 combinations (2^3) of treating each mediator as a confounder:
#   Model A: none (total effect — correct under DAG)
#   Model B: + MDA only
#   Model C: + Activity only
#   Model D: + Duration only
#   Model E: + Activity + MDA
#   Model F: + MDA + Duration
#   Model G: + Activity + Duration
#   Model H: + Activity + MDA + Duration (most conservative — all as confounders)
#
# Note: Duration sits downstream of Activity in the DAG (Activity -> Duration ->
# Infection), so Models D/F/G/H additionally test whether this sub-pathway
# structure matters.
# =============================================================================
message("\n=== SENSITIVITY 4: DAG STRUCTURAL UNCERTAINTY (ALL MEDIATORS) ===\n")

# Model A: total effect (already fitted as model_travel_inf / model_travel_tot)
# formula: infected_bin ~ travel_frequency + age_class + Sex_bin + Location

# Model B: + MDA only
model_s4_mda <- load_or_fit(
  "model_s4_mda.rds",
  formula = infected_bin ~ travel_frequency + MDA + age_class + Sex_bin + Location,
  data = data, family = bernoulli(), prior = prior_bern
)

# Model C: + Activity only (previously model_travel_act_adjusted)
model_s4_act <- load_or_fit(
  "model_s4_act.rds",
  formula = infected_bin ~ travel_frequency + ActNEW + age_class + Sex_bin + Location,
  data = data, family = bernoulli(), prior = prior_bern
)

# Model D: + Duration only
model_s4_dur <- load_or_fit(
  "model_s4_dur.rds",
  formula = infected_bin ~ travel_frequency + DurMin + age_class + Sex_bin + Location,
  data = data, family = bernoulli(), prior = prior_bern
)

# Model E: + Activity + MDA
model_s4_act_mda <- load_or_fit(
  "model_s4_act_mda.rds",
  formula = infected_bin ~ travel_frequency + ActNEW + MDA + age_class + Sex_bin + Location,
  data = data, family = bernoulli(), prior = prior_bern
)

# Model F: + MDA + Duration
model_s4_mda_dur <- load_or_fit(
  "model_s4_mda_dur.rds",
  formula = infected_bin ~ travel_frequency + MDA + DurMin + age_class + Sex_bin + Location,
  data = data, family = bernoulli(), prior = prior_bern
)

# Model G: + Activity + Duration
model_s4_act_dur <- load_or_fit(
  "model_s4_act_dur.rds",
  formula = infected_bin ~ travel_frequency + ActNEW + DurMin + age_class + Sex_bin + Location,
  data = data, family = bernoulli(), prior = prior_bern
)

# Model H: + Activity + MDA + Duration (most conservative)
model_s4_all <- load_or_fit(
  "model_s4_all.rds",
  formula = infected_bin ~ travel_frequency + ActNEW + MDA + DurMin + age_class + Sex_bin + Location,
  data = data, family = bernoulli(), prior = prior_bern
)

# ---- Extract daily travel OR (92x vs never) for each model ------------------
extract_daily_or <- function(model, model_label, mediators_conditioned) {
  d    <- as_draws_df(model)
  vals <- d[[daily_coef]]
  tibble(
    model                = model_label,
    mediators            = mediators_conditioned,
    log_or_med           = median(vals),
    log_or_lo            = as.numeric(quantile(vals, 0.025)),
    log_or_hi            = as.numeric(quantile(vals, 0.975)),
    or_med               = exp(median(vals)),
    or_lo                = exp(as.numeric(quantile(vals, 0.025))),
    or_hi                = exp(as.numeric(quantile(vals, 0.975))),
    cri_clears_null      = exp(as.numeric(quantile(vals, 0.025))) > 1
  )
}

dag_sensitivity <- bind_rows(
  extract_daily_or(model_travel_tot, "A: Total effect (correct DAG)",        "None"),
  extract_daily_or(model_s4_mda,     "B: + MDA",                             "MDA"),
  extract_daily_or(model_s4_act,     "C: + Activity",                        "Activity"),
  extract_daily_or(model_s4_dur,     "D: + Duration",                        "Duration"),
  extract_daily_or(model_s4_act_mda, "E: + Activity + MDA",                  "Activity + MDA"),
  extract_daily_or(model_s4_mda_dur, "F: + MDA + Duration",                  "MDA + Duration"),
  extract_daily_or(model_s4_act_dur, "G: + Activity + Duration",             "Activity + Duration"),
  extract_daily_or(model_s4_all,     "H: + Activity + MDA + Duration",       "All mediators")
)

# Order for display (Model A at top)
dag_sensitivity$model <- factor(dag_sensitivity$model, levels = rev(dag_sensitivity$model))

cat("=== Travel effect (daily vs never): across all structural assumptions ===\n")
print(dag_sensitivity %>%
        select(model, mediators, or_med, or_lo, or_hi, log_or_med, cri_clears_null),
      n = Inf)

# ---- Forest plot (daily travel OR across all 8 models) ----------------------
p_dag <- ggplot(dag_sensitivity,
                aes(x = or_med, xmin = or_lo, xmax = or_hi, y = model,
                    colour = cri_clears_null)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_errorbarh(height = 0.25, linewidth = 0.8) +
  geom_point(size = 3) +
  scale_colour_manual(values = c("TRUE" = "#2166ac", "FALSE" = "#d6604d"),
                      labels = c("TRUE" = "CrI above null", "FALSE" = "CrI crosses null"),
                      name   = NULL) +
  scale_x_log10(breaks = c(0.5, 1, 2, 4, 8)) +
  labs(x     = "Odds ratio for daily travel vs. never (log scale, 95% CrI)",
       y     = NULL,
       title = "Sensitivity 4: DAG structural uncertainty",
       subtitle = "Effect of treating each mediator as a confounder — daily travel vs. never") +
  theme_minimal(base_size = 13) +
  theme(plot.background  = element_rect(fill = "white", colour = NA),
        legend.position  = "bottom",
        panel.grid.minor = element_blank())

ggsave("Figures/dag_structural_sensitivity.png", p_dag, width = 9, height = 6, dpi = 300)
cat("Figure saved: Figures/dag_structural_sensitivity.png\n")

write.csv(dag_sensitivity %>% mutate(model = as.character(model)),
          "Figures/dag_structural_sensitivity.csv", row.names = FALSE)


# =============================================================================
# COMBINED SUMMARY TABLE
# =============================================================================
cat("\n\n======================================================\n")
cat("SUMMARY OF SENSITIVITY ANALYSES\n")
cat("======================================================\n\n")

cat("--- SENSITIVITY 1: E-values (all key associations) ---\n")
print(evalue_supp %>%
        select(outcome, exposure, effect_type, effect_est, cri_lo, cri_hi,
               evalue_estimate, evalue_inner_cri), n = Inf)

cat("\n--- SENSITIVITY 2a: Measurement error (DurMin) ---\n")
cat("Direction of bias: TOWARD NULL (observed effect is conservative)\n")
print(me_sensitivity %>%
        select(reliability_ratio, pct_measurement_error, or_per_min_corrected,
               log_or_corrected_lo, log_or_corrected_hi), n = Inf)

cat("\n--- SENSITIVITY 3: Spillover / substitution bounding ---\n")
print(spillover_combined %>%
        select(type, spillover_frac, diff_med_pct, diff_lo_pct, diff_hi_pct), n = Inf)

cat("\n--- SENSITIVITY 4: DAG structural uncertainty (daily travel OR, all 8 models) ---\n")
print(dag_sensitivity %>%
        mutate(model = as.character(model)) %>%
        select(model, mediators, or_med, or_lo, or_hi, cri_clears_null), n = Inf)

cat("\nAll output tables saved to Figures/ directory.\n")
cat("All models saved to", model_dir, "/ directory for reloading.\n")

# =============================================================================
# MCMC CONVERGENCE DIAGNOSTICS TABLE
# Reports R-hat and effective sample sizes (bulk + tail) for all key
# parameters across the six main E-value models and the duration model.
# Criteria (following Vehtari et al. 2021): R-hat < 1.01; ESS > 400.
# =============================================================================
message("\n=== MCMC CONVERGENCE DIAGNOSTICS ===\n")

extract_diagnostics <- function(model, model_label) {
  s <- posterior::summarise_draws(
    posterior::as_draws(model),
    mean, sd,
    ~posterior::quantile2(.x, probs = c(0.025, 0.975)),
    posterior::default_convergence_measures()
  )
  # Keep only population-level (b_) and intercept parameters
  s <- s[grepl("^b_|^Intercept", s$variable), ]
  s$model <- model_label
  s[, c("model", "variable", "mean", "sd", "q2.5", "q97.5", "rhat",
        "ess_bulk", "ess_tail")]
}

mcmc_diag <- bind_rows(
  extract_diagnostics(model_travel_inf,  "Travel → Infection"),
  extract_diagnostics(model_act_inf,     "Activity → Infection"),
  extract_diagnostics(model_mda_inf,     "MDA → Infection"),
  extract_diagnostics(model_travel_burd, "Travel → Burden"),
  extract_diagnostics(model_act_burd,    "Activity → Burden"),
  extract_diagnostics(model_mda_burd,    "MDA → Burden"),
  extract_diagnostics(model_dur_tot,     "Duration → Infection")
)

mcmc_diag <- mcmc_diag %>%
  mutate(across(c(mean, sd, q2.5, q97.5), ~round(.x, 3)),
         rhat     = round(rhat, 4),
         ess_bulk = round(ess_bulk),
         ess_tail = round(ess_tail),
         flag     = ifelse(rhat >= 1.01 | ess_bulk < 400 | ess_tail < 400,
                           "CHECK", "ok"))

cat("=== MCMC convergence diagnostics (R-hat < 1.01; ESS > 400 required) ===\n")
print(mcmc_diag, n = Inf)

n_flags <- sum(mcmc_diag$flag == "CHECK", na.rm = TRUE)
if (n_flags == 0) {
  cat("\nAll parameters pass convergence criteria (R-hat < 1.01, ESS > 400).\n")
} else {
  cat("\nWARNING:", n_flags, "parameter(s) flagged — inspect trace plots.\n")
  print(mcmc_diag %>% filter(flag == "CHECK"), n = Inf)
}

write.csv(mcmc_diag, "Figures/mcmc_convergence_diagnostics.csv", row.names = FALSE)
cat("Convergence table saved: Figures/mcmc_convergence_diagnostics.csv\n")
