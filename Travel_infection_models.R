#Full model of travel frequncy and infection

rm(list = ls())
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/HPGH/_PhD/Jess_data/data_manipulation")

library(dagitty)
library(ggdag)
library(brms)
library(dplyr)
library(ggplot2)

#install data
#
data <-read.csv("Travel_prevEDIT.csv")
nrow(data) #586

#select the columns I want

# data <- data %>%
#   select(Location, age_class, Age, Sex_bin, mean_sm, mean_sm_plus, infection,infected_bin, Lived, ActNEW, CL, DurMin, MDA,travel_frequency, travel_bin)

nrow(data)

#How many NAs
#
colSums(is.na(data)) #none yaay

data$ActNEW <- as.factor(data$ActNEW)
data$ActNEW <- relevel(data$ActNEW, ref = "None") 
data$Actlong <- as.factor(data$Actlong)
data$Actlong <- relevel(data$Actlong, ref = "None") 
data$Testing_loc <- as.factor(data$Testing_loc)
data$Testing_loc <- relevel(data$Testing_loc, ref = "NAAWA") 
data$infection <- as.factor(data$infection)
data$infection <- relevel(data$infection, ref = "0") 

#Exposure freq outcome infection
#Direct = age, sex, location, MDA, Activity
#Total = age, sex, location

#TOTAL
#Infection YN
priors <- get_custom_prior("bernoulli")

#What about travel as a categorical variable?

data$travel_frequency <- as.factor(data$travel_frequency)

#TOTAL
#Do infection yes or no first 
Freq_inf_bin <- brm(infected_bin~ travel_frequency+age_class+Sex_bin+Location,
                    data=data,
                    family="bernoulli",
                    prior=priors,
                    iter = 10000, warmup = 3000, chains = 4, cores = 8,
                    control = list(max_treedepth = 15))

#model fit
pp_check(Freq_inf_bin)  # shows dens_overlay plot by default
pp_check(Freq_inf_bin, type = "error_hist", ndraws = 11)
pp_check(Freq_inf_bin, type = "stat_2d") 

summary(Freq_inf_bin) #

plot(Freq_inf_bin)

#DIRECT
#Do infection yes or no first 
Freq_inf_bin_dir <- brm(infected_bin~ travel_frequency+age_class+Sex_bin+Location+ActNEW+MDA,
                    data=data,
                    family="bernoulli",
                    prior=priors,
                    iter = 10000, warmup = 3000, chains = 4, cores = 8,
                    control = list(max_treedepth = 25))

#model fit
pp_check(Freq_inf_bin_dir)  # shows dens_overlay plot by default
pp_check(Freq_inf_bin_dir, type = "error_hist", ndraws = 11)
pp_check(Freq_inf_bin, type = "stat_2d") 

summary(Freq_inf_bin_dir) #

#then mean_sm only TOTAL
filt <- data %>%
  filter(mean_sm >0)#%>%
  #mutate(mean_sm= round(mean_sm))
nrow(filt) #64

priors <- get_custom_prior("gamma")

Freq_inf_mean <- brm(mean_sm~ travel_frequency+age_class+Sex_bin+Location,
                    data=filt,
                    family="gamma",
                    prior = priors,
                    iter = 10000, warmup = 3000, chains = 4, cores = 8,
                    control = list(max_treedepth = 15))

#model fit
pp_check(Freq_inf_mean )  # shows dens_overlay plot by default
pp_check(Freq_inf_mean , type = "error_hist", ndraws = 11)
pp_check(Freq_inf_mean , type = "stat_2d") 

summary(Freq_inf_mean ) #no effect



#DIRCT
Freq_inf_mean_dir <- brm(mean_sm~ travel_frequency+age_class+Sex_bin+Location+MDA+ActNEW,
                     data=filt,
                     family="gamma",
                     prior=priors,
                     iter = 10000, warmup = 3000, chains = 4, cores = 8,
                     control = list(max_treedepth = 15))

#model fit
pp_check(Freq_inf_mean_dir )  # shows dens_overlay plot by default
pp_check(Freq_inf_mean_dir , type = "error_hist", ndraws = 11)
pp_check(Freq_inf_mean_dir , type = "stat_2d") 

summary(Freq_inf_mean_dir ) #no effect


#make plots

library(tidybayes)
library(ggplot2)
library(dplyr)
library(tidyr)
library(viridis)


#what are the names in my model
variables(Freq_inf_bin)


draws <- Freq_inf_bin %>%
  spread_draws(b_travel_frequency1, b_travel_frequency2, b_travel_frequency3,
               b_travel_frequency6, b_travel_frequency13, b_travel_frequency26,
               b_travel_frequency92)

#reshape for plotting
draws_long <- draws %>%
  pivot_longer(cols = starts_with("b_travel_frequency"),
               names_to = "travel_frequency",
               values_to = "estimate") %>%
  mutate(travel_frequency = case_when(
    travel_frequency == "b_travel_frequency1" ~ "1",
    travel_frequency == "b_travel_frequency2" ~ "2",
    travel_frequency == "b_travel_frequency3" ~ "3",
    travel_frequency == "b_travel_frequency6" ~ "6",
    travel_frequency == "b_travel_frequency13" ~ "13",
    travel_frequency == "b_travel_frequency26" ~ "26",
    travel_frequency == "b_travel_frequency92" ~ "92",
    TRUE ~ travel_frequency
  ))

# Make travel_frequency a factor with the correct order
draws_long <- draws_long %>%
  mutate(travel_frequency = factor(travel_frequency,
                                   levels = c("92", "26", "13", "6", "3", "2", "1")))


# Now plot posterior distributions
ggplot(draws_long, aes(x = estimate, fill = travel_frequency)) +
  stat_halfeye(slab_alpha = 0.6, .width = c(0.66, 0.95)) +
  facet_wrap(~ travel_frequency, scales = "fixed", ncol=1) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(
    x = "Effect on log-odds of infection",
    y = NULL,
    fill = "Days travelled to the lake\nin the last 3 months",  # <- legend title
    title = NULL
  ) +
  theme_minimal() +
  scale_fill_viridis_d(option = "viridis")+
  theme(
    axis.text.x = element_text(size = 14, face="bold"),   # x-axis text
    axis.text.y = element_blank(),                        # y-axis text
                     # x-axis title
    axis.title.x = element_text(size = 14, face="bold"),                         # y-axis title
    strip.text = element_blank(),                           # facet titles (village names)
    plot.title = element_text(size = 18, hjust = 0.5, face="bold"),              # main plot title
    legend.text = element_text(size = 14),                          # legend labels
    legend.title = element_text(size = 14, face="bold"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank()
   
  ) +
  xlim(-3,3)

ggsave("Figures/post_dist_freq_inf_tot.png", width = 10, height = 10, dpi = 300)


# Now plot trace
# 
library(bayesplot)

table(data$travel_frequency)

# Extract the posterior draws as a matrix
posterior_matrix <- as.matrix(Freq_inf_bin$fit)

# Now plot trace plots for travel frequency terms only
# Make parameter a factor with correct levels
posterior_long <- posterior_long %>%
  mutate(parameter = factor(parameter, levels = c(
    "b_travel_frequency92",
    "b_travel_frequency26",
    "b_travel_frequency13",
    "b_travel_frequency6",
    "b_travel_frequency3",
    "b_travel_frequency2",
    "b_travel_frequency1"
  )))




library(posterior)
library(dplyr)
library(tidyr)
library(ggplot2)

# Step 1: Convert to array format
posterior_array <- as_draws_array(Freq_inf_bin)

# Step 2: Subset just the travel_frequency parameters
param_names <- names(travel_colours)
posterior_array_sub <- posterior_array[, , param_names, drop = FALSE]

# Step 3: Convert to data frame with proper chain/iteration info
posterior_long <- posterior_array_sub %>%
  posterior::as_draws_df() %>%
  mutate(.chain = rep(1:4, each = 7000)) %>%         # 4 chains x 7000 iterations post-warmup
  mutate(.iteration = rep(1:7000, times = 4)) %>%
  pivot_longer(cols = all_of(param_names),
               names_to = "parameter",
               values_to = "value")


ggplot(posterior_long, aes(x = .iteration, y = value, colour = parameter)) +
  geom_line(alpha = 0.7) +
  scale_color_viridis_d(option = "viridis") +
  facet_wrap(~ parameter, ncol = 1, scales = "free_y") +
  labs(
    x = "Iteration",
    y = NULL,
    colour = "Days travelled",
    title = NULL
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    strip.text = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    text = element_text(size = 14, face="bold"),  
    axis.text.x = element_text(face="bold")# x-axis text
  )

ggsave("Figures/trace_freq_inf_tot.png", width = 4, height = 12, dpi = 300)

#why are people who travel once a week less liekly to be infected?

filter <- data%>%
  filter(travel_frequency =="13")

table(filter$Actlong)

names(filter)


######
levels_order <- c("1", "2", "3", "6", "13", "26", "92")  # ascending order

#Total effect on infection ######
plot_posterior_distributions(
  model = Freq_inf_bin,
  param_prefix = "travel_frequency",
  levels_order = levels_order,
  label_title = "Days travelled to the lake\nin the last 3 months",
  output_file = "Figures/test/post_dist_freq_inf_tot.png",
  xlab = "Effect on log-odds of infection"
)

plot_trace(
  model = Freq_inf_bin,
  param_prefix = "travel_frequency",
  levels_order = levels_order,
  output_file = "Figures/test/trace_freq_inf_tot.png"
)

#Direct effect on infection #######
plot_posterior_distributions(
  model = Freq_inf_bin_dir,
  param_prefix = "travel_frequency",
  levels_order = levels_order,
  label_title = "Days travelled to the lake\nin the last 3 months",
  output_file = "Figures/test/post_dist_freq_inf_dir.png",
  xlab = "Effect on log-odds of infection"
)

plot_trace(
  model = Freq_inf_bin_dir,
  param_prefix = "travel_frequency",
  levels_order = levels_order,
  output_file = "Figures/test/trace_freq_inf_dir.png"
)


#Total effect on EPG #######
plot_posterior_distributions(
  model = Freq_inf_mean,
  param_prefix = "travel_frequency",
  levels_order = levels_order,
  label_title = "Days travelled to the lake\nin the last 3 months",
  output_file = "Figures/test/post_dist_freq_epg_tot.png",
  xlab = "Effect on log mean eggs per gram"
)

plot_trace(
  model = Freq_inf_mean,
  param_prefix = "travel_frequency",
  levels_order = levels_order,
  output_file = "Figures/test/trace_freq_epg_tot.png"
)

#Direcct effect on EPG #######
plot_posterior_distributions(
  model = Freq_inf_mean_dir,
  param_prefix = "travel_frequency",
  levels_order = levels_order,
  label_title = "Days travelled to the lake\nin the last 3 months",
  output_file = "Figures/test/post_dist_freq_epg_dir.png",
  xlab = "Effect on log mean eggs per gram"
)

plot_trace(
  model = Freq_inf_mean_dir,
  param_prefix = "travel_frequency",
  levels_order = levels_order,
  output_file = "Figures/test/trace_freq_epg_dir.png"
)
