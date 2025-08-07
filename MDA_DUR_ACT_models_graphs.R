#models and graph for MDA, duration and activity effects on Infection prevalence and burden



rm(list = ls())
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/HPGH/_PhD/Jess_data/data_manipulation")

library(brms)
library(tidybayes)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggdist)
library(viridis)

#install data
#
data <-read.csv("Travel_prevEDIT.csv")
nrow(data) #585


colSums(is.na(data)) #none yaay

data$ActNEW <- as.factor(data$ActNEW)
data$ActNEW <- relevel(data$ActNEW, ref = "None") 
data$Actlong <- as.factor(data$Actlong)
data$Actlong <- relevel(data$Actlong, ref = "None") 
data$Location<- as.factor(data$Location)
data$Location <- relevel(data$Location, ref = "Naawa") 

#travel as a categorical variable?

data$travel_frequency <- as.factor(data$travel_frequency)

#Load function called priors.R

#============================================================================#
###MDA ###
#============================================================================#
#TOTAL and DIRECT
#Age, location, Travel freq

#Prevalence ####

MDA <- brm(infected_bin ~ MDA + age_class+Location+travel_frequency,
           data=data,
           family="bernoulli",
           prior=c(prior(normal(0, 10), class = "b"),
           prior(normal(-1.5, 2.5), class = "Intercept")),
           iter = 10000, warmup = 3000, chains = 4, cores = 8,
           control = list(max_treedepth = 15))


get_prior(MDA)

#model fit
pp_check(MDA)  # shows dens_overlay plot by default
pp_check(MDA, type = "error_hist", ndraws = 11)
pp_check(MDA, type = "stat_2d") 

summary(MDA) #

#Intensity ####

filt <- data%>%
  filter(mean_sm >0)

MDA_mean_sm <- brm(mean_sm ~ MDA + age_class+Location+travel_frequency,
                   data=filt,
                   family="gamma",
                  prior=  c(
                    prior(normal(2.5, 1), class = "Intercept"),
                    prior(normal(0, 10), class = "b"),
                    prior(exponential(1), class = "shape")),
                   iter = 10000, warmup = 3000, chains = 4, cores = 8,
                   control = list(max_treedepth = 15))


get_prior(MDA_mean_sm)

#model fit
pp_check(MDA_mean_sm )  # shows dens_overlay plot by default
pp_check(MDA_mean_sm , type = "error_hist", ndraws = 11)
pp_check(MDA_mean_sm , type = "stat_2d") 

summary(MDA_mean_sm ) #

#============================================================================#
###DURATION ###
#============================================================================#


#TOTAL and DIRECT
#Age, Sex, Location, Activity 

#Prevalence #####


Dur <- brm(infected_bin ~ DurMin+age_class+Location+Sex_bin+ActNEW,
           data=data,
           family="bernoulli",
           prior=c(prior(normal(0, 10), class = "b"),
                   prior(normal(-1.5, 2.5), class = "Intercept")),
           iter = 10000, warmup = 3000, chains = 4, cores = 8,
           control = list(max_treedepth = 15))

#model fit
pp_check(Dur)  # shows dens_overlay plot by default
pp_check(Dur, type = "error_hist", ndraws = 11)
pp_check(Dur, type = "stat_2d") 

summary(Dur) #

#
#Intensity #####

filt <- data%>%
  filter(mean_sm >0)


Dur_mean_sm <- brm(mean_sm ~ DurMin + age_class+Location+Sex_bin+ActNEW,
                   data=filt,
                   family="gamma",  
                   prior=  c(
                     prior(normal(2.5, 1), class = "Intercept"),
                     prior(normal(0, 10), class = "b"),
                     prior(exponential(1), class = "shape")),
                   iter = 10000, warmup = 3000, chains = 4, cores = 8,
                   control = list(max_treedepth = 15))

#model fit
pp_check(Dur_mean_sm )  # shows dens_overlay plot by default
pp_check(Dur_mean_sm , type = "error_hist", ndraws = 11)
pp_check(Dur_mean_sm , type = "stat_2d") 

summary(Dur_mean_sm ) #

#what about intensity with only people who spend time in the lake

filt <- data%>%
  filter(mean_sm >0)

filt2 <- filt%>%
  filter(DurMin >0)

Dur_mean_sm2 <- brm(mean_sm ~ DurMin + age_class+Location+Sex_bin+ActNEW,
                   data=filt2,
                   family="gamma",
                   prior=  c(
                     prior(normal(2.5, 1), class = "Intercept"),
                     prior(normal(0, 10), class = "b"),
                     prior(exponential(1), class = "shape")),
                   iter = 10000, warmup = 3000, chains = 4, cores = 8,
                   control = list(max_treedepth = 15))

#model fit
pp_check(Dur_mean_sm2 )  # shows dens_overlay plot by default
pp_check(Dur_mean_sm2 , type = "error_hist", ndraws = 11)
pp_check(Dur_mean_sm2 , type = "stat_2d") 

summary(Dur_mean_sm2 ) #



#============================================================================#
###Activity ###
#============================================================================#
#
#Prevalence####

#TOTAL 
#Age, Sex, Location, Travel freq


Act <- brm(infected_bin ~ ActNEW+age_class+Location+Sex_bin+travel_frequency,
           data=data,
           family="bernoulli",
           prior=c(prior(normal(0, 10), class = "b"),
                   prior(normal(-1.5, 2.5), class = "Intercept")),
           iter = 10000, warmup = 3000, chains = 4, cores = 8,
           control = list(max_treedepth = 15))

#model fit
pp_check(Act)  # shows dens_overlay plot by default
pp_check(Act, type = "error_hist", ndraws = 11)
pp_check(Act, type = "stat_2d") 

summary(Act) #



#DIRECT

Act_Dur <- brm(infected_bin ~ ActNEW+age_class+Location+Sex_bin+travel_frequency+DurMin,
               data=data,
               family="bernoulli",
               prior=c(prior(normal(0, 10), class = "b"),
                       prior(normal(-1.5, 2.5), class = "Intercept")),
               iter = 10000, warmup = 3000, chains = 4, cores = 8,
               control = list(max_treedepth = 15))


summary(Act_Dur) #



#Intesnity ####

#TOTAL
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

#model fit
pp_check(Act_mean_sm )  # shows dens_overlay plot by default
pp_check(Act_mean_sm , type = "error_hist", ndraws = 11)
pp_check(Act_mean_sm , type = "stat_2d") 

summary(Act_mean_sm ) #



#DIRECT
Act_mean_sm_dur <- brm(mean_sm ~ ActNEW+ age_class+Location+Sex_bin+travel_frequency+DurMin,
                   data=filt,
                   family="gamma",
                   prior=  c(
                     prior(normal(2.5, 1), class = "Intercept"),
                     prior(normal(0, 10), class = "b"),
                     prior(exponential(1), class = "shape")),
                   iter = 10000, warmup = 3000, chains = 4, cores = 8,
                   control = list(max_treedepth = 15))

get_prior(Act_mean_sm_dur)


#model fit
pp_check(Act_mean_sm_dur )  # shows dens_overlay plot by default
pp_check(Act_mean_sm_dur , type = "error_hist", ndraws = 11)
pp_check(Act_mean_sm_dur , type = "stat_2d") 

summary(Act_mean_sm_dur ) #

#=========================================================================================================================
#Create posterior distribution and trace plots for Total effects of MDA, Duration and Activity on infection prevalence ####
##=========================================================================================================================
#Which models?
# MDA
# Dur
# Act
# 
# Extract posterior draws ####
library(tidybayes)
library(tidyverse)

# 1. Extract posterior draws

# MDA (from model `MDA`)
mda_draws <- MDA %>%
  spread_draws(b_MDA) %>%
  mutate(parameter = "MDA") %>%
  rename(.value = b_MDA)

# Duration (from model `Dur`)
duration_draws <- Dur %>%
  spread_draws(b_DurMin) %>%
  mutate(parameter = "Duration") %>%
  rename(.value = b_DurMin)

# Activity (from model `Act`)
activity_draws <- Act %>%
  spread_draws(b_ActNEWDomestic,
               b_ActNEWOccupation,
               b_ActNEWRecreation,
               b_ActNEWTradeORVisit) %>%
  pivot_longer(
    cols = starts_with("b_ActNEW"),
    names_to = "parameter",
    values_to = ".value"
  ) %>%
  mutate(parameter = recode(parameter,
                            "b_ActNEWDomestic" = "Activity: Domestic",
                            "b_ActNEWOccupation" = "Activity: Occupation",
                            "b_ActNEWRecreation" = "Activity: Recreation",
                            "b_ActNEWTradeORVisit" = "Activity: Trade/Visit"
  ))

# 2. Combine all draws
all_draws <- bind_rows(mda_draws,  activity_draws)

#Set the order for the plots

# Define the order from top (first) to bottom (last)
param_levels <- c(
  "Activity: Domestic",
  "Activity: Occupation",
  "Activity: Recreation",
  "Activity: Trade/Visit",
  "MDA"
)

# Apply it to the full dataset
all_draws$parameter <- factor(all_draws$parameter, levels = param_levels)


# 3. Assign unique .draw IDs per parameter (required for proper density)
all_draws <- all_draws %>%
  group_by(parameter) %>%
  mutate(.draw = row_number()) %>%
  ungroup()

# 4. Plot posterior distributions
posterior_plot <- ggplot(all_draws, aes(x = .value, y = parameter, fill = parameter)) +
  stat_halfeye(
    .width = c(0.66, 0.95),
    slab_alpha = 0.6,
    point_interval = "mean_qi"
  ) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_fill_viridis_d(option = "plasma") +
  theme_minimal() +
  labs(
    x = "Estimated effect (log-odds of infection)",
    y = NULL,
    fill = NULL
  ) +
  theme(
    text = element_text(size = 20, face="bold"),
    legend.position = "none",
    
    plot.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA)
  )



# 5. Display plot
posterior_plot

#Separet plot for duration
duration_draws$y <- 1  # assign numeric y-position

duration_plot <- ggplot(duration_draws, aes(x = .value, y = y)) +
  stat_halfeye(
    aes(fill = "Duration"),  # to still fill by parameter
    .width = c(0.66, 0.95),
    slab_alpha = 0.6,
    point_interval = "mean_qi",
    height = 0.2
  ) +
  scale_y_continuous(
    breaks = 1,
    labels = "Duration",
    expand = c(0.1, 0.1)  # controls spacing above/below
  ) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c("Duration" = "#666666")) +
  labs(
    x = "Effect of duration on log-odds of infection",
    y = NULL,
    fill = NULL
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 20, face="bold"),
    legend.position = "none",
    strip.text = element_blank(),
    panel.spacing = unit(0, "lines"),
    plot.margin = margin(0, 0, 0, 0),
    
    plot.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA)
  )




duration_plot

ggsave("Figures/posterior_plot_weak_prior.png", posterior_plot, width = 8, height = 12, dpi = 300)
ggsave("Figures/duration_plot_weak_prior.png", duration_plot, width = 5, height = 2, dpi = 300)


# Set the desired top-to-bottom order for facets
param_levels <- c(
  "MDA",
  "Activity: Trade/Visit",
  "Activity: Recreation",
  "Activity: Occupation",
  "Activity: Domestic"
)

# Apply the factor level order
all_draws$parameter <- factor(all_draws$parameter, levels = param_levels)

# Plot
trace_plot <- ggplot(all_draws, aes(x = .iteration, y = .value, colour = parameter)) +
  geom_line(alpha = 0.6) +
  facet_wrap(~parameter, ncol = 1, scales = "free_y", strip.position = "left", as.table = FALSE) +
  scale_color_viridis_d(option = "plasma") +
  labs(x = "Iteration", y = NULL, colour = NULL) +
  theme_minimal() +
  theme(
    text = element_text(size = 20, face="bold"),
    axis.text.y = element_blank(),
    strip.text = element_blank(),
    panel.grid.major.y = element_blank(),
    
    plot.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA)
  )


trace_plot#Combine plots

ggsave("Figures/trace_plot_posts_weak_prior.png", trace_plot, width = 8, height = 12, dpi = 300)

#trace plot for duration
#
duration_trace <- ggplot(duration_draws, aes(x = .iteration, y = .value)) +
  geom_line(alpha = 0.6, colour = "#666666") +
  labs(x = "Iteration", y = NULL) +
  theme_minimal() +
  theme(
    text = element_text(size = 20, face="bold"),
    axis.text.y = element_blank(),
    strip.text = element_blank(),
    panel.grid.major.y = element_blank(),
    
    plot.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.margin = margin(0, 10, 0, 10)  # tighten margins
  ) +
  ylim(min(duration_draws$.value), max(duration_draws$.value))

duration_trace

ggsave("Figures/duration_trace_weak_prior.png", duration_trace, height = 2, width = 4, dpi = 300)



#=========================================================================================================================
#Create posterior distribution and trace plots for Total effects of MDA, Duration and Activity on infection INTENSITY ####
##=========================================================================================================================
#Intensity plots
#
## 1. Extract posterior draws for egg count models

# MDA (egg count)
mda_sm_draws <- MDA_mean_sm %>%
  spread_draws(b_MDA) %>%
  mutate(parameter = "MDA") %>%
  rename(.value = b_MDA)

# Duration (egg count)
duration_sm_draws <- Dur_mean_sm %>%
  spread_draws(b_DurMin) %>%
  mutate(parameter = "Duration") %>%
  rename(.value = b_DurMin)

# Activity (egg count)
activity_sm_draws <- Act_mean_sm %>%
  spread_draws(b_ActNEWDomestic,
               b_ActNEWOccupation,
               b_ActNEWRecreation,
               b_ActNEWTradeORVisit) %>%
  pivot_longer(
    cols = starts_with("b_ActNEW"),
    names_to = "parameter",
    values_to = ".value"
  ) %>%
  mutate(parameter = recode(parameter,
                            "b_ActNEWDomestic" = "Activity: Domestic",
                            "b_ActNEWOccupation" = "Activity: Occupation",
                            "b_ActNEWRecreation" = "Activity: Recreation",
                            "b_ActNEWTradeORVisit" = "Activity: Trade/Visit"
  ))

# 2. Combine all draws (excluding Duration for separate plot)
all_draws_sm <- bind_rows(mda_sm_draws, activity_sm_draws)

# Set factor levels in reverse for consistent top-down order in both plots
param_levels_sm <- c(
  "Activity: Domestic",
  "Activity: Occupation",
  "Activity: Recreation",
  "Activity: Trade/Visit",
  "MDA"
)

param_levels_sm_rev <- rev(param_levels_sm)

# Apply reversed order globally
all_draws_sm$parameter <- factor(all_draws_sm$parameter, levels = param_levels_sm_rev)


# 4. Assign unique draw IDs
all_draws_sm <- all_draws_sm %>%
  group_by(parameter) %>%
  mutate(.draw = row_number()) %>%
  ungroup()

# 5. Posterior plot
posterior_plot_sm <- ggplot(all_draws_sm, aes(x = .value, y = parameter, fill = parameter)) +
  stat_halfeye(.width = c(0.66, 0.95), slab_alpha = 0.6, point_interval = "mean_qi") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_fill_viridis_d(option = "plasma") +
  labs(x = "Estimated effect (log mean egg count)", y = NULL, fill = NULL) +
  theme_minimal() +
  theme(
    text = element_text(size = 20, face="bold"),
    legend.position = "none",
    plot.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA)
  )

posterior_plot_sm
ggsave("Figures/posterior_plot_eggcount2_weak_prior.png", posterior_plot_sm, width = 8, height = 12, dpi = 300)

# 6. Duration posterior for egg count
duration_sm_draws$y <- 1

duration_plot_sm <- ggplot(duration_sm_draws, aes(x = .value, y = y)) +
  stat_halfeye(
    aes(fill = "Duration"),
    .width = c(0.66, 0.95),
    slab_alpha = 0.6,
    point_interval = "mean_qi",
    height = 0.2
  ) +
  scale_y_continuous(breaks = 1, labels = "Duration", expand = c(0.1, 0.1)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c("Duration" = "#666666")) +
  labs(x = "Estimated effect (log mean egg count)", y = NULL, fill = NULL) +
  theme_minimal() +
  theme(
    text = element_text(size = 20, face="bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    strip.text = element_blank(),
    panel.spacing = unit(0, "lines"),
    plot.margin = margin(0, 0, 0, 0),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA)
  )

duration_plot_sm

ggsave("Figures/duration_plot_eggcount2_weak_prior.png", duration_plot_sm, width = 4, height = 2, dpi = 300)

# 7. Trace plot for activity + MDA
param_levels_sm_rev <- rev(param_levels_sm)
all_draws_sm$parameter <- factor(all_draws_sm$parameter, levels = param_levels_sm_rev)

trace_plot_sm <- ggplot(all_draws_sm, aes(x = .iteration, y = .value, colour = parameter)) +
  geom_line(alpha = 0.6) +
  facet_wrap(~parameter, ncol = 1, scales = "free_y", strip.position = "left", as.table = FALSE) +
  scale_color_viridis_d(option = "plasma") +
  labs(x = "Iteration", y = NULL, colour = NULL) +
  theme_minimal() +
  theme(
    text = element_text(size = 14),
    axis.text.y = element_blank(),
    strip.text = element_blank(),
    panel.grid.major.y = element_blank(),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA)
  )
trace_plot_sm
ggsave("Figures/trace_plot_eggcount2.png", trace_plot_sm, width = 8, height = 12, dpi = 300)

# 8. Trace plot for Duration
duration_trace_sm <- ggplot(duration_sm_draws, aes(x = .iteration, y = .value)) +
  geom_line(alpha = 0.6, colour = "#666666") +
  labs(x = "Iteration", y = NULL) +
  theme_minimal() +
  theme(
    text = element_text(size = 14),
    axis.text.y = element_blank(),
    strip.text = element_blank(),
    panel.grid.major.y = element_blank(),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.margin = margin(0, 10, 0, 10)
  ) +
  ylim(min(duration_sm_draws$.value), max(duration_sm_draws$.value))

duration_trace_sm

ggsave("Figures/duration_trace_eggcount2_weak_prior.png", duration_trace_sm, height = 2, width = 4, dpi = 300)


#=========================================================================================================================
#Supplementary graphs for direct effect of activity so including duration ####
##=========================================================================================================================

#pick the corresponding colours

library(viridisLite)
library(scales)

# Get the "magma" palette with N colours
magma_colours <- viridisLite::plasma(n = 5)

# Show them
scales::show_col(magma_colours)

activity_colours <- c(
  "Activity: Domestic"     = "#F0F921FF",
  "Activity: Occupation"   = "#F89441FF",
  "Activity: Recreation"   = "#CC4678FF",
  "Activity: Trade/Visit"  = "#7E03A8FF"
)


#Prevalence #######
#
## Extract draws from prevalence model
act_prev_draws <- Act_Dur %>%
  spread_draws(b_ActNEWDomestic, b_ActNEWOccupation, b_ActNEWRecreation, b_ActNEWTradeORVisit) %>%
  pivot_longer(
    cols = starts_with("b_ActNEW"),
    names_to = "parameter",
    values_to = ".value"
  ) %>%
  mutate(parameter = recode(parameter,
                            "b_ActNEWDomestic" = "Activity: Domestic",
                            "b_ActNEWOccupation" = "Activity: Occupation",
                            "b_ActNEWRecreation" = "Activity: Recreation",
                            "b_ActNEWTradeORVisit" = "Activity: Trade/Visit"
  ))

# Factor order for plotting
param_levels <- c("Activity: Domestic", "Activity: Occupation", "Activity: Recreation", "Activity: Trade/Visit")
act_prev_draws$parameter <- factor(act_prev_draws$parameter, levels = rev(param_levels))

# Posterior plot (prevalence)
posterior_prev <- ggplot(act_prev_draws, aes(x = .value, y = parameter, fill = parameter)) +
  stat_halfeye(.width = c(0.66, 0.95), slab_alpha = 0.6, point_interval = "mean_qi") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = activity_colours) +
  labs(x = "Effect on log-odds of infection", y = NULL, fill = NULL) +
  theme_minimal() +  
  theme(
    text = element_text(size = 14),
    legend.position = "none",
    plot.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA)
  )

posterior_prev 
ggsave("Figures/SUP_posterior_plot_activity_prevalence.png", posterior_prev, width = 8, height = 12, dpi = 300)

# Trace plot (prevalence)
trace_prev <- ggplot(act_prev_draws, aes(x = .iteration, y = .value, colour = parameter)) +
  geom_line(alpha = 0.6) +
  facet_wrap(~parameter, ncol = 1, scales = "free_y", as.table = FALSE) +
  scale_color_manual(values = activity_colours) +
  labs(x = "Iteration", y = NULL, colour = NULL) +
  theme_minimal() +
  theme(
    text = element_text(size = 14),
    axis.text.y = element_blank(),
    strip.text = element_blank(),
    panel.grid.major.y = element_blank(),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA)
  )


trace_prev
ggsave("Figures/SUP_trace_plot_activity_prevalence.png", trace_prev, width = 8, height = 12, dpi = 300)

#Intensity ########
#
# Extract draws from egg burden model
act_burden_draws <- Act_mean_sm_dur %>%
  spread_draws(b_ActNEWDomestic, b_ActNEWOccupation, b_ActNEWRecreation, b_ActNEWTradeORVisit) %>%
  pivot_longer(
    cols = starts_with("b_ActNEW"),
    names_to = "parameter",
    values_to = ".value"
  ) %>%
  mutate(parameter = recode(parameter,
                            "b_ActNEWDomestic" = "Activity: Domestic",
                            "b_ActNEWOccupation" = "Activity: Occupation",
                            "b_ActNEWRecreation" = "Activity: Recreation",
                            "b_ActNEWTradeORVisit" = "Activity: Trade/Visit"
  ))

act_burden_draws$parameter <- factor(act_burden_draws$parameter, levels = rev(param_levels))

# Posterior plot (egg burden)
posterior_burden <- ggplot(act_burden_draws, aes(x = .value, y = parameter, fill = parameter)) +
  stat_halfeye(.width = c(0.66, 0.95), slab_alpha = 0.6, point_interval = "mean_qi") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = activity_colours) +
  labs(x = "Effect on log mean egg count", y = NULL, fill = NULL) +
  theme_minimal() +
  theme(
    text = element_text(size = 14),
    legend.position = "none",
    plot.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA)
  )

posterior_burden
ggsave("Figures/SUP_posterior_plot_activity_burden.png", posterior_burden, width = 8, height = 12, dpi = 300)

# Trace plot (egg burden)
trace_burden <- ggplot(act_burden_draws, aes(x = .iteration, y = .value, colour = parameter)) +
  geom_line(alpha = 0.6) +
  facet_wrap(~parameter, ncol = 1, scales = "free_y", as.table = FALSE) +
  scale_color_manual(values = activity_colours) +
  labs(x = "Iteration", y = NULL, colour = NULL) +
  theme_minimal() +
  theme(
    text = element_text(size = 14),
    axis.text.y = element_blank(),
    strip.text = element_blank(),
    panel.grid.major.y = element_blank(),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA)
  )

trace_burden
ggsave("Figures/SUP_trace_plot_activity_burden.png", trace_burden, width = 8, height = 12, dpi = 300)

