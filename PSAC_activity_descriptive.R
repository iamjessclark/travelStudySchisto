# Descriptive analysis: infected PSAC by activity type
# Generated in response to reviewer comment asking for data on
# number of young children infected and their activity type.
# This is a descriptive summary only — cell sizes are too small
# for formal modelling of activity effects within PSAC.

library(dplyr)

data <- read.csv("Travel_prevEDIT.csv")

# --- PSAC: n total, n infected, % infected by activity type ---
psac_by_activity <- data %>%
  filter(age_class == "PSAC") %>%
  group_by(ActNEW) %>%
  summarise(
    n_total    = n(),
    n_infected = sum(infected_bin == 1, na.rm = TRUE),
    pct_infected = round(100 * n_infected / n_total, 1),
    .groups = "drop"
  ) %>%
  arrange(desc(n_total))

print(psac_by_activity)

# --- Of all infected PSAC: breakdown by activity type ---
psac_infected_activity <- data %>%
  filter(age_class == "PSAC", infected_bin == 1) %>%
  count(ActNEW, name = "n_infected") %>%
  mutate(pct_of_infected_psac = round(100 * n_infected / sum(n_infected), 1)) %>%
  arrange(desc(n_infected))

print(psac_infected_activity)

# Summary stats
cat("Total PSAC:", sum(data$age_class == "PSAC"), "\n")
cat("Total infected PSAC:", sum(data$age_class == "PSAC" & data$infected_bin == 1), "\n")
cat("Overall PSAC prevalence:",
    round(100 * mean(data$infected_bin[data$age_class == "PSAC"]), 1), "%\n")
