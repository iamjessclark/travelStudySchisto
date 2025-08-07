#create maps of travel site points

# Load libraries
library(ggmap)
library(ggplot2)
library(dplyr)
library(ggspatial)
library(ggrepel)


library(MetBrewer)
MetBrewer::display_all(colorblind_only = T, n=12) 

library(paletteer)
paletteer_d("MetBrewer::Signac")


#####low risk

# Define locations with specific coordinates (these are the home coordinates found on the travel survey csv and converted to decimal using chat gpt)
locations <- data.frame(
  Villages = c("B", "Nw", "Nam", "S", "Nan"),
  lattitude = c(0.3399, 0.3496, 0.1792, 0.3722, 0.3761),
  longitude = c(33.1605, 33.1914, 32.9369, 33.2101, 33.2442)
)


# Compute the center of all locations
center_lat <- mean(locations$lattitude)
center_lon <- mean(locations$longitude)

# Register Google API key
ggmap::register_google(key="")

# Adjust the center of the map based on mean latitude and longitude
uganda_map <- get_map(location = c(lon = center_lon, lat = center_lat), zoom = 11, maptype = "terrain")

# Plot the map with flagged locations and custom colors for each place
ggmap(uganda_map) +
  geom_point(data = locations, aes(x = longitude, y = lattitude, color = Villages), size = 5) +
  scale_color_manual(values = c(
    "B" = "#663171FF",
    "Nw" = "#CF3A36FF",
    "Nam" = "#EA7428FF",
    "S" = "#E2998AFF",
    "Nan" = "#0C7156FF"
  )) +
  theme_minimal() +
  labs(
    x = "Longitude",
    y = "Latitude",
    title = "Village Locations")+
  theme(
    axis.text.x = element_text(hjust = 1, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 14, face="bold"),
    axis.title.y = element_text(size = 14, face="bold"),
    strip.text = element_text(size = 14, face="bold"),
    plot.title = element_text(size = 18, hjust = 0.5, face="bold"),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14, face="bold"),
    plot.background = element_rect(fill = "white", colour = NA)
  )

ggsave("Figures/village_locations.png", width = 12, height = 8, dpi = 300)

#Zoom out to show Uganda

library(ggmap)
# Adjust the center of the map
uganda_map <- get_map(location = "Uganda", zoom = 7, maptype = "terrain")

# Plot the map with flagged locations and custom colors for each place
ggmap(uganda_map) +
  geom_point(data = locations, aes(x = longitude, y = lattitude, colour=Villages), size = 2) +
  #geom_label_repel(data = locations, aes(x = longitude, y = lattitude, label = place),
  #box.padding = 0.35, point.padding = 0.5, segment.color = 'grey50') +
  scale_color_manual(values = c(
    "B" = "#663171FF",
    "Nw" = "#CF3A36FF",
    "Nam" = "#EA7428FF",
    "S" = "#E2998AFF",
    "Nan" = "#0C7156FF"
  )) +
  theme_minimal() +
  labs(
    x = "Longitude",
    y = "Latitude",
    title = "Village Locations")+
  theme(
    axis.text.x = element_text(hjust = 1, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 14, face="bold"),
    axis.title.y = element_text(size = 14, face="bold"),
    strip.text = element_text(size = 14, face="bold"),
    plot.title = element_text(size = 18, hjust = 0.5, face="bold"),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14, face="bold"),
    plot.background = element_rect(fill = "white", colour = NA)
  )


#####high risk check if it is Namukuma or Ssese that we are doing for low risk and change accordingly

# Define locations with specific coordinates
locationsHIGH <- data.frame(
  High_Risk_Locations = c("Bugoba", "Buwagajjo", "Kalega", "Nanso B"),
  lattitude = c(0.3177, 0.334028, 0.347611, 0.368528),
  longitude = c(33.1754, 33.197028, 33.223361, 33.259)
)

# Register Google API key
ggmap::register_google(key="")

# Adjust the center of the map
uganda_map <- get_map(location = c(lon = 33.1, lat = 0.2), zoom = 11, maptype = "terrain")

# Plot the map with flagged locations and custom colors for each place
ggmap(uganda_map) +
  geom_point(data = locationsHIGH, aes(x = longitude, y = lattitude, color = High_Risk_Locations), size = 3.5) +
  #geom_label_repel(data = locations, aes(x = longitude, y = lattitude, label = place),
  #box.padding = 0.35, point.padding = 0.5, segment.color = 'grey50') +
  scale_color_manual(values = c("Bugoba" = "pink", "Buwagajjo" = "dark Blue", 
                                "Kalega" = "red", "Nanso B" = "yellow")) +
  theme_minimal() +
  ylim(0.2,0.4)+
  xlim(33.0,33.3)


###zoom in high risk


# Adjust the center of the map
uganda_map <- get_map(location = c(lon = 33.2, lat = 0.2), zoom = 11, maptype = "terrain")



# Plot the map with flagged locations and custom colors for each place
ggmap(uganda_map) +
  geom_point(data = locationsHIGH, aes(x = longitude, y = lattitude, color = High_Risk_Locations), size = 3.5) +
  #geom_label_repel(data = locations, aes(x = longitude, y = lattitude, label = place),
  #box.padding = 0.35, point.padding = 0.5, segment.color = 'grey50') +
  scale_color_manual(values = c("Bugoba" = "pink", "Buwagajjo" = "dark Blue", 
                                "Kalega" = "red", "Nanso B" = "yellow")) +
  theme_minimal() +
  ylim(0.3,0.4)+
  xlim(33.1,33.3)

#Zoom in both high and low

# Low risk data frame
locations <- data.frame(
  Location = c("Banga II", "Naawa", "Sesse", "Nanso A"),
  Risk = "Low Risk",
  latitude = c(0.3399, 0.3496, 0.3722, 0.3761),
  longitude = c(33.1605, 33.1914, 33.2101, 33.2442)
)

# High risk data frame
locationsHIGH <- data.frame(
  Location = c("Bugoba", "Buwagajjo", "Kalega", "Nanso B"),
  Risk = "High Risk",
  latitude = c(0.3177, 0.334028, 0.347611, 0.368528),
  longitude = c(33.1754, 33.197028, 33.223361, 33.259)
)

# Combine both data frames
all_locations <- rbind(
  data.frame(Location = locationsHIGH$Location, Risk = locationsHIGH$Risk, latitude = locationsHIGH$latitude, longitude = locationsHIGH$longitude)
,   data.frame(Location = locations$Location, Risk = locations$Risk, latitude = locations$latitude, longitude = locations$longitude)
)


# Create the map
ggmap(uganda_map) +
  geom_point(data = all_locations, aes(x = longitude, y = latitude, color = Location), size = 3.5) +
  scale_color_manual(values = c(
    "Bugoba" = "pink", "Buwagajjo" = "dark Blue", "Kalega" = "red", "Nanso B" = "yellow",
    "Banga II" = "purple", "Naawa" = "Blue", "Sesse" = "dark red", "Nanso A" = "orange")) +
  theme_minimal() +
  ylim(0.3, 0.4) + xlim(33.1, 33.3) +
  guides(color = guide_legend(title = "Regions", override.aes = list(size = 4))) +
  facet_wrap(~Risk, ncol = 1)+
  theme(text = element_text(face="bold", size = 20))+
  ylab("Lattitude")+
  xlab("Longitude")

#plot all datapoints on the same map

# Plot the map with flagged locations and custom colors for each place
ggmap(uganda_map) +
  geom_point(data = all_locations, aes(x = longitude, y = latitude, color = Location), size = 3.5) +
  #geom_label_repel(data = locations, aes(x = longitude, y = lattitude, label = place),
  #box.padding = 0.35, point.padding = 0.5, segment.color = 'grey50') +
  scale_color_manual(values = c("Bugoba" = "pink", "Buwagajjo" = "dark Blue", 
                                "Kalega" = "red", "Nanso B" = "yellow",  
                                "Banga II" = "purple", "Naawa" = "Blue", 
                                "Sesse" = "dark red", "Nanso A" = "orange")) +
  theme_minimal() +
  ylim(0.3,0.4)+
  xlim(33.1,33.3)


