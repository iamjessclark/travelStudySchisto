#maps of where particpants siad their "home water contact" site was

setwd("~/Library/Mobile Documents/com~apple~CloudDocs/HPGH/_PhD/Jess_data/data_manipulation")

library(dplyr)
library(ggmap)
library(measurements)

# Read in the data
locations_raw <- read.csv("locations.csv", stringsAsFactors = FALSE)

# Clean names
dms_to_decimal <- function(dms) {
  # Remove whitespace and split into components
  dms <- trimws(dms)
  matches <- regmatches(dms, regexec("([0-9]+)°([0-9]+)'([0-9.]+)\"?([NSEW])", dms))
  
  if (length(matches[[1]]) < 5) return(NA)  # Return NA if format doesn't match
  
  deg <- as.numeric(matches[[1]][2])
  min <- as.numeric(matches[[1]][3])
  sec <- as.numeric(matches[[1]][4])
  dir <- matches[[1]][5]
  
  dec <- deg + min / 60 + sec / 3600
  
  # Flip sign if West or South
  if (dir %in% c("S", "W")) dec <- -dec
  
  return(dec)
}

locations_clean <- locations_raw %>%
  mutate(
    name = trimws(gsub("Water.*|stream.*|landing.*", "", name, ignore.case = TRUE)),
    latitude = sapply(latitude_dms, dms_to_decimal),
    longitude = sapply(longitude_dms, dms_to_decimal)
  )


center_lat <- mean(locations_clean$latitude, na.rm = TRUE)
center_lon <- mean(locations_clean$longitude, na.rm = TRUE)

# Register Google API key
ggmap::register_google(key="")

# Compute bounding box
bbox <- make_bbox(lon = longitude, lat = latitude, data = locations_clean, f = 0.1)  # f = 0.1 adds a bit of padding

# Download map using bounding box
uganda_map <- get_map(location = bbox, maptype = "terrain", source = "google")

# Plot all points without cutting any off
ggmap(uganda_map) +
  geom_point(data = locations_clean, aes(x = longitude, y = latitude, color = name), alpha = 0.6) +
  theme_minimal() +
  labs(title = "Mapped Survey Coordinates") +
  guides(color = guide_legend(title = "Water Contact"))




