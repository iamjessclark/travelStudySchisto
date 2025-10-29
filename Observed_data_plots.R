# =========================
# Real-data plots (no models)
# =========================

# Packages
library(tidyverse)
library(binom)      # for Wilson CIs
library(scales)     # for nice axis formatting
library(patchwork)  # to combine plots (optional but handy)
library(purrr)
library(viridis)

#code generated using real data and AI to write code. All checked against real data for mistakes.

#-------------------------
# 0) Helper functions
#-------------------------

# Prevalence + 95% CI for a categorical exposure
prev_by_group <- function(df, group_var, infected_var = "infected_bin") {
  df %>%
    mutate(.group = forcats::fct_infreq(.data[[group_var]])) %>%
    group_by(.group) %>%
    summarise(
      n = n(),
      y = sum(.data[[infected_var]] == 1, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      ci_list = map2(y, n, ~ binom::binom.wilson(.x, .y)),  # list of data.frames
      lwr = map_dbl(ci_list, ~ .x$lower),
      upr = map_dbl(ci_list, ~ .x$upper),
      prev = y / n
    ) %>%
    select(.group, n, prev, lwr, upr)
}


prev_by_bins <- function(df, x, infected_var = "infected_bin", n_bins = 12, breaks = NULL) {
  x_sym <- ensym(x)
  
  # pull non-missing values for break construction
  v <- df %>% pull(!!x_sym)
  v <- v[!is.na(v)]
  
  # construct breaks
  if (is.null(breaks)) {
    qs <- quantile(v, probs = seq(0, 1, length.out = n_bins + 1), na.rm = TRUE)
    br <- unique(as.numeric(qs))
    # fallback if quantiles collapse (few unique values)
    if (length(br) < 3) {
      br <- pretty(v, n = min(n_bins, max(1, length(unique(v)) - 1)))
    }
  } else {
    br <- sort(unique(breaks))
  }
  
  # if still not enough breaks, return empty tibble gracefully
  if (length(br) < 2) {
    return(tibble(
      bin = integer(), xmid = numeric(), prev = numeric(),
      lwr = numeric(), upr = numeric(), n = integer(),
      xmin = numeric(), xmax = numeric()
    ))
  }
  
  df %>%
    filter(!is.na(!!x_sym)) %>%
    mutate(
      bin = cut(!!x_sym, breaks = br, include.lowest = TRUE, right = TRUE, labels = FALSE)
    ) %>%
    group_by(bin) %>%
    summarise(
      n = n(),
      y = sum(.data[[infected_var]] == 1, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      xmin = br[pmax(1, bin)],
      xmax = br[pmin(length(br) - 1, bin) + 1L],
      xmid = (xmin + xmax) / 2,
      ci_df = map2(y, n, ~ binom::binom.wilson(.x, .y)),
      lwr = map_dbl(ci_df, ~ .x$lower),
      upr = map_dbl(ci_df, ~ .x$upper),
      prev = y / n
    ) %>%
    select(bin, xmid, prev, lwr, upr, n, xmin, xmax)
}
# A minimal, consistent theme
theme_pub <- function(base_size = 12) {
  theme_classic(base_size = base_size) +
    theme(
      panel.border = element_blank(),
      axis.title = element_text(face = "bold"),
      strip.background = element_rect(fill = "grey95", colour = NA),
      strip.text = element_text(face = "bold")
    )
}

#-------------------------
# 1) Tidy/ensure types
#-------------------------
data_plot <- data %>%
  mutate(
    travel_frequency = as.factor(travel_frequency),
    MDA = as.factor(MDA),         
    ActNEW = as.factor(ActNEW),
    infected_bin = as.integer(infected_bin),
    mean_sm = as.numeric(mean_sm),
    DurMin = as.numeric(DurMin)
  )

#-------------------------
# 2) Infection (binary) vs categorical exposures
#    A) travel_frequency
#-------------------------
prev_freq <- prev_by_group(data_plot, "travel_frequency")

# Define the correct travel-frequency order
freq_levels <- c("0", "1", "2", "3", "6", "13", "26", "92")

# Make sure both datasets use that order
prev_freq$.group <- factor(prev_freq$.group, levels = freq_levels)
data_plot$travel_frequency <- factor(data_plot$travel_frequency, levels = freq_levels)

p_inf_freq <- ggplot(prev_freq, aes(x = .group, y = prev, fill = .group)) +
  geom_col(colour = "grey40", alpha = 0.8) +
  geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.15) +
  geom_jitter(
    data = data_plot,
    aes(x = travel_frequency, y = infected_bin),
    width = 0.15, height = 0.03, alpha = 0.25, size = 0.6,
    inherit.aes = FALSE
  ) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
  scale_fill_viridis_d(option = "viridis") +
  labs(x = "Travel frequency (visits in 3 months)", y = "Infection prevalence (observed)") +
  theme_pub() +
  theme(legend.position = "none")
p_inf_freq

#-------------------------
#    B) MDA
#-------------------------
prev_mda <- prev_by_group(data_plot, "MDA")

# (Optional) rename for clarity
prev_mda <- prev_by_group(data_plot, "MDA") |>
  dplyr::rename(MDA = .group)

# make sure MDA is a factor (0 = no, 1 = yes)
data_plot$MDA <- factor(data_plot$MDA, levels = c("0", "1"), labels = c("No", "Yes"))
prev_mda$MDA <- factor(prev_mda$MDA, levels = c("0", "1"), labels = c("No", "Yes"))

p_inf_MDA <- ggplot(prev_mda, aes(x = MDA, y = prev, fill = MDA)) +
  geom_col(colour = "grey40", alpha = 0.8) +
  geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.15) +
  geom_jitter(
    data = data_plot,
    aes(x = MDA, y = infected_bin),
    width = 0.15, height = 0.03, alpha = 0.25, size = 0.6,
    inherit.aes = FALSE
  ) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
  scale_fill_viridis_d(option = "plasma") +   # 🔥 plasma palette
  labs(
    x = "MDA participation",
    y = "Infection prevalence (observed)"
  ) +
  theme_pub() +
  theme(legend.position = "none")

#-------------------------
#    C) Activity
#-------------------------
prev_act <- prev_by_group(data_plot, "ActNEW")

# Ensure factor order (edit to your preferred display order)
act_levels <- c("None","Domestic","Recreation","TradeORVisit","Occupation")
data_plot$ActNEW <- factor(data_plot$ActNEW, levels = act_levels)

# Prevalence + Wilson CIs
prev_act <- prev_by_group(data_plot, "ActNEW") |>
  dplyr::rename(ActNEW = .group)
prev_act$ActNEW <- factor(prev_act$ActNEW, levels = act_levels)

# Plot
p_inf_act <- ggplot(prev_act, aes(x = ActNEW, y = prev, fill = ActNEW)) +
  geom_col(colour = "grey40", alpha = 0.8) +
  geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.15) +
  geom_jitter(
    data = data_plot,
    aes(x = ActNEW, y = infected_bin),
    width = 0.15, height = 0.03, alpha = 0.25, size = 0.6,
    inherit.aes = FALSE
  ) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
  scale_fill_viridis_d(option = "plasma") +  # viridis palette for activity
  labs(x = "Activity type", y = "Infection prevalence (observed)") +
  theme_pub() +
  theme(legend.position = "none")


#-------------------------
# 3) Infection (binary) vs duration (continuous)
#    Binned prevalence + rugs of raw DurMin + jittered 0/1
#-------------------------
prev_dur <- prev_by_bins(data_plot, DurMin, breaks = seq(0, 720, by = 60))

p_inf_dur <- ggplot(prev_dur, aes(x = xmid, y = prev)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.15, fill = "grey70") +
  geom_line(size = 0.8) +
  geom_point(size = 1.2) +
  geom_rug(
    data = data_plot,
    aes(x = DurMin),
    sides = "b", alpha = 0.25,
    inherit.aes = FALSE
  ) +
  geom_jitter(
    data = data_plot,
    aes(x = DurMin, y = infected_bin),
    height = 0.03, width = 0, alpha = 0.15, size = 0.5,
    inherit.aes = FALSE
  ) +
  labs(x = "Minutes in Lake Victoria (DurMin)", y = "Infection prevalence (binned, observed)") +
  theme_pub()

p_inf_dur

# If DurMin is very skewed, consider uncommenting:
#p_inf_dur <- p_inf_dur + scale_x_continuous(trans = "log1p")
# Choose behaviourally-meaningful bins (edit as you like)
breaks <- c(0, 1, 10, 30, 60, Inf)
labels <- c("0", "1–10", "11–30", "31–60", ">60")

data_binned <- data %>%
  mutate(
    DurBin = cut(DurMin, breaks = breaks, include.lowest = TRUE, right = TRUE, labels = labels)
  ) %>%
  filter(!is.na(DurBin))

prev_bin <- data_binned %>%
  group_by(DurBin) %>%
  summarise(
    n = n(),
    y = sum(infected_bin == 1, na.rm = TRUE),
    prev = y / n,
    lwr = binom::binom.wilson(y, n)$lower,
    upr = binom::binom.wilson(y, n)$upper,
    .groups = "drop"
  )


p_inf_dur2 <- ggplot(prev_bin, aes(x = DurBin, y = prev)) +
  geom_col(fill = "grey85") +
  geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.15) +
  geom_jitter(data = data_binned,
              aes(x = DurBin, y = infected_bin),
              width = 0.15, height = 0.03, alpha = 0.20, size = 0.6) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0,1)) +
  labs(x = "Minutes in Lake Victoria (binned)", y = "Infection prevalence (observed)") +
  theme_classic(base_size = 12)

ggplot(prev_bin, aes(x = DurBin, y = prev)) +
  geom_pointrange(aes(ymin = lwr, ymax = upr), size = 0.4) +
  # optional: show raw 0/1s jittered within bins
  geom_jitter(data = data_binned,
              aes(x = DurBin, y = infected_bin),
              width = 0.15, height = 0.03, alpha = 0.20, size = 0.6) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0,1)) +
  labs(x = "Minutes in Lake Victoria (binned)", y = "Infection prevalence (observed)") +
  theme_classic(base_size = 12)

#-------------------------
# 4) Intensity (EPG) vs categorical exposures
#    Use half-violin + box (or just boxplot) with raw points
#-------------------------

# Common EPG transform for skew/zeros
data_plot <- data_plot %>% mutate(EPG_log1p = log1p(mean_sm))

# Ensure correct factor order for x-axis
freq_levels <- c("0", "1", "2", "3", "6", "13", "26", "92")
data_plot$travel_frequency <- factor(data_plot$travel_frequency, levels = freq_levels)

# Log-transform EPG to reduce skew
data_plot <- data_plot %>%
  mutate(EPG_log1p = log1p(EPG))

# Violin + box + jitter plot
p_epg_freq <- ggplot(data_plot, aes(x = travel_frequency, y = EPG_log1p, fill = travel_frequency)) +
  geom_violin(trim = FALSE, colour = "grey70", alpha = 0.7) +
  geom_jitter(width = 0.15, height = 0, alpha = 0.25, size = 0.6) +
  scale_fill_viridis_d(option = "viridis") +   # keep viridis palette for travel frequency
  labs(
    x = "Travel frequency (visits in 3 months)",
    y = "Infection intensity (EPG, log1p)"
  ) +
  theme_pub() +
  theme(legend.position = "none")

p_epg_freq
p_epg_mda <- ggplot(data_plot, aes(x = MDA, y = EPG_log1p)) +
  geom_violin(fill = "grey90", colour = "grey70", scale = "count") +
  geom_boxplot(width = 0.15, outlier.shape = NA) +
  geom_jitter(width = 0.15, height = 0, alpha = 0.25, size = 0.6) +
  scale_y_continuous(name = "EPG (log1p scale)", breaks = pretty_breaks()) +
  labs(x = "MDA participation") +
  theme_pub()

p_epg_act <- ggplot(data_plot, aes(x = ActNEW, y = EPG_log1p, fill = ActNEW)) +
  geom_violin(trim = FALSE, colour = "grey70", alpha = 0.7, scale = "count") +

  geom_jitter(width = 0.15, height = 0, alpha = 0.25, size = 0.6) +
  scale_fill_viridis_d(option = "plasma", direction = -1) +  # same order
  labs(x = "Activity type", y = "EPG (log1p)") +
  theme_pub() +
  theme(legend.position = "none")
p_epg_act
#-------------------------
# 5) Intensity (EPG) vs duration (continuous)
#    Scatter + smooth; add rugs for DurMin & EPG_log1p
#-------------------------
p_epg_dur <- ggplot(data_plot, aes(x = DurMin, y = EPG_log1p)) +
  geom_point(alpha = 0.25, size = 0.7) +
  geom_smooth(method = "loess", se = TRUE, span = 0.9) +
  geom_rug(sides = "b", alpha = 0.25) +
  labs(x = "Minutes in Lake Victoria (DurMin)", y = "EPG (log1p scale)") +
  theme_pub()
# Optionally: + scale_x_continuous(trans = "log1p")

#-------------------------
# 6) (Optional) assemble figures
#-------------------------
fig_inf_cats <- p_inf_freq + p_inf_mda + p_inf_act + plot_layout(ncol = 3)
fig_inf_cont <- p_inf_dur2
fig_epg_cats <- p_epg_freq + p_epg_mda + p_epg_act + plot_layout(ncol = 3)
fig_epg_cont <- p_epg_dur

# Print to screen
fig_inf_cats
fig_inf_cont
fig_epg_cats
fig_epg_cont

#-------------------------
# 7) Save (optional)
#-------------------------
ggsave("fig_infection_vs_categorical.png", fig_inf_cats, width = 12, height = 4, dpi = 300)
ggsave("fig_infection_vs_duration.png",   fig_inf_cont, width = 6, height = 4, dpi = 300)
ggsave("fig_intensity_vs_categorical.png", fig_epg_cats, width = 12, height = 4, dpi = 300)
ggsave("fig_intensity_vs_duration.png",    fig_epg_cont, width = 6, height = 4, dpi = 300)
