# Descriptive analysis: activity type by age class and infection status
# Generated in response to reviewer comment asking for data on
# number of children and adults infected and their activity type.
# This is a descriptive summary only — cell sizes may be too small
# for formal modelling of activity effects within age classes.

library(dplyr)

data <- read.csv("Travel_prevEDIT.csv")

# --- Combined table: all age classes, activity type, infection ---
# pct_of_infected is the % of all infected within that age class
activity_table <- data %>%
  group_by(age_class, ActNEW) %>%
  summarise(
    n_total      = n(),
    n_infected   = sum(infected_bin == 1, na.rm = TRUE),
    pct_infected = round(100 * n_infected / n_total, 1),
    .groups = "drop"
  ) %>%
  group_by(age_class) %>%
  mutate(
    pct_of_infected = round(100 * n_infected / sum(n_infected), 1)
  ) %>%
  ungroup() %>%
  arrange(age_class, desc(n_total))

print(activity_table)

write.csv(activity_table, "activity_table.csv")
# Summary stats by age class
data %>%
  group_by(age_class) %>%
  summarise(
    n_total    = n(),
    n_infected = sum(infected_bin == 1, na.rm = TRUE),
    prevalence = round(100 * n_infected / n_total, 1),
    .groups = "drop"
  ) %>%
  print()

# --- Validation: do "None" activity participants = non-travellers? ---
# Cross-tabulate ActNEW == "None" against travel_bin to check consistency.
# Expectation: None activity should align with travel_bin == 0 (did not travel).
cat("\n--- Validation: ActNEW 'None' vs travel_bin (all participants) ---\n")
print(table(ActNEW = data$ActNEW, travel_bin = data$travel_bin))

cat("\n--- Validation: ActNEW 'None' vs travel_bin (PSAC only) ---\n")
print(table(
  ActNEW     = data$ActNEW[data$age_class == "PSAC"],
  travel_bin = data$travel_bin[data$age_class == "PSAC"]
))

# The small number with None + travel_bin==1 are local residents (same location)
# who visit the site rarely but reported no specific activity purpose.
cat("\nNone-activity participants with travel_bin==1: check residence\n")
none_travelers <- data[
  !is.na(data$ActNEW) & data$ActNEW == "None" & data$travel_bin == 1,
]
cat("n =", nrow(none_travelers), "\n")
cat("Resident at same location:\n")
print(table(none_travelers$Is.that.the.same.location.as.here., useNA = "always"))
cat("Visit frequency:\n")
print(table(none_travelers$Freq, useNA = "always"))
cat("Age class:\n")
print(table(none_travelers$age_class))
