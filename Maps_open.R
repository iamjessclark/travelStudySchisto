# Maps using pure OpenStreetMap tiles via maptiles (no API key required)
# Replaces Google Maps API used in Maps.R to comply with copyright requirements.

#
# Changes from Maps.R:
#   - ggmap/Google replaced with maptiles + tidyterra + sf
#   - No API key or account registration needed
#   - ggmap() replaced with ggplot() + geom_spatraster_rgb()
#   - All point/colour/theme/ggsave code unchanged


library(ggplot2)
library(dplyr)
library(ggspatial)
library(ggrepel)
library(MetBrewer)
library(paletteer)
library(maptiles)   # OSM tile fetching
library(tidyterra)  # geom_spatraster_rgb for plotting tiles
library(sf)         # bbox handling

MetBrewer::display_all(colorblind_only = T, n = 12)
paletteer_d("MetBrewer::Signac")

# Helper: build an sf bbox and fetch OSM tiles
osm_tile <- function(left, bottom, right, top, zoom = 11) {
  bbox_sf <- st_bbox(
    c(xmin = left, ymin = bottom, xmax = right, ymax = top),
    crs = st_crs(4326)
  ) %>% st_as_sfc()
  get_tiles(bbox_sf, provider = "OpenStreetMap", zoom = zoom, crop = TRUE)
}

# =============================================================================
# LOW RISK VILLAGES
# =============================================================================

locations <- data.frame(
  Villages  = c("B", "Nw", "Nam", "S", "Nan"),
  lattitude = c(0.3399, 0.3496, 0.1792, 0.3722, 0.3761),
  longitude = c(33.1605, 33.1914, 32.9369, 33.2101, 33.2442)
)

center_lat <- mean(locations$lattitude)
center_lon <- mean(locations$longitude)

village_tile <- osm_tile(
  left   = min(locations$longitude) - 0.08,
  bottom = min(locations$lattitude)  - 0.20,  # extended south to show Lake Victoria
  right  = max(locations$longitude)  + 0.08,
  top    = max(locations$lattitude)  + 0.05,
  zoom   = 10
)

ggplot() +
  geom_spatraster_rgb(data = village_tile) +
  geom_point(data = locations, aes(x = longitude, y = lattitude, color = Villages), size = 5) +
  scale_color_manual(values = c(
    "B"   = "#663171FF",
    "Nw"  = "#CF3A36FF",
    "Nam" = "#EA7428FF",
    "S"   = "#E2998AFF",
    "Nan" = "#0C7156FF"
  )) +
  theme_minimal() +
  labs(x = "Longitude", y = "Latitude") +
  theme(
    axis.text.x  = element_text(hjust = 1, size = 14),
    axis.text.y  = element_text(size = 14),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    strip.text   = element_text(size = 14, face = "bold"),
    plot.title   = element_text(size = 18, hjust = 0.5, face = "bold"),
    legend.text  = element_text(size = 14),
    legend.title = element_text(size = 14, face = "bold"),
    plot.background = element_rect(fill = "white", colour = NA)
  )


ggsave("Figures/village_locations.png", width = 12, height = 8, dpi = 300) 

# =============================================================================
# UGANDA OVERVIEW
# =============================================================================

uganda_tile <- osm_tile(
  left = 29.5, bottom = -1.5, right = 35.1, top = 4.2,
  zoom = 7
)

ggplot() +
  geom_spatraster_rgb(data = uganda_tile) +
  geom_point(data = locations, aes(x = longitude, y = lattitude, colour = Villages), size = 2) +
  scale_color_manual(values = c(
    "B"   = "#663171FF",
    "Nw"  = "#CF3A36FF",
    "Nam" = "#EA7428FF",
    "S"   = "#E2998AFF",
    "Nan" = "#0C7156FF"
  )) +
  theme_minimal() +
  labs(x = "Longitude", y = "Latitude") +
  theme(
    axis.text.x  = element_text(hjust = 1, size = 14),
    axis.text.y  = element_text(size = 14),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    strip.text   = element_text(size = 14, face = "bold"),
    plot.title   = element_text(size = 18, hjust = 0.5, face = "bold"),
    legend.text  = element_text(size = 14),
    legend.title = element_text(size = 14, face = "bold"),
    plot.background = element_rect(fill = "white", colour = NA)
  )

ggsave("Figures/Uganda.png", width = 12, height = 8, dpi = 300) 
