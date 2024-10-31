# Set working directory
wd_str <- file.path("F:", "replication_repository_YJUEC_D_24_00004R2")

install_if_not_exists <- function(package) {
  if (!require(package, character.only = TRUE)) {
    install.packages(package, dependencies = TRUE)
    library(package, character.only = TRUE)
  }
}
install_if_not_exists("data.table")
install_if_not_exists("ggplot2")
install_if_not_exists("ggmap")
install_if_not_exists("sf")
install_if_not_exists("RColorBrewer")
install_if_not_exists("svglite")
install_if_not_exists("dplyr")

setwd(wd_str)

# Read and filter data
dt <- fread(file.path("processed_data", "df_map.csv"))
dt <- dt[F_bar %in% c(0.5, 0.6, 0.9, 1.25, 2.0), ]
dt[, F_bar := as.factor(F_bar)]

# Convert to sf object
points_sf <- st_as_sf(dt, coords = c("XCoord", "YCoord"), crs = 2263)

# Check if geometries are valid before transformation
if (any(!st_is_valid(points_sf))) {
  points_sf <- st_make_valid(points_sf)
}

# Transform coordinates
points_sf <- st_transform(points_sf, 4326)

# Load NYC boroughs shapefile
nyc_boroughs <- st_read(file.path("additional_files", "borough_boundaries",
                                  "geo_export.shp"))

# Load each borough's shapefile and check for valid geometries
load_and_check <- function(filepath) {
  shp <- st_read(filepath)
  if (any(!st_is_valid(shp))) {
    shp <- st_make_valid(shp)
  }
  return(shp)
}

bk_map <- load_and_check(file.path("raw_data", "MAPPLUTO",
                                   "citywide", "BKMapPLUTO.shp"))
bx_map <- load_and_check(file.path("raw_data", "MAPPLUTO",
                                   "citywide", "BXMapPLUTO.shp"))
mn_map <- load_and_check(file.path("raw_data", "MAPPLUTO",
                                   "citywide", "MNMapPLUTO.shp"))
qn_map <- load_and_check(file.path("raw_data", "MAPPLUTO",
                                   "citywide", "QNMapPLUTO.shp"))
si_map <- load_and_check(file.path("raw_data", "MAPPLUTO",
                                   "citywide", "SIMapPLUTO.shp"))

nyc_map <- rbind(bk_map, bx_map, mn_map, qn_map, si_map)

# Check if the combined map has valid geometries
if (any(!st_is_valid(nyc_map))) {
  nyc_map <- st_make_valid(nyc_map)
}

color_palette <- brewer.pal(n = length(unique(dt$F_bar)), name = "Set1")

my_plot <- ggplot() +
  geom_sf(data = nyc_boroughs, fill = "lightgray", color = "black") +
  geom_sf(data = points_sf, aes(color = as.factor(F_bar)),
          size = 1.8, alpha = 0.5) +
  scale_color_manual(values = color_palette) +  # Customize the colors as needed
  theme_minimal() +
  labs(title = "Figure 8: Location of Sample Properties by FAR Limit",
       color = "FAR Limit")

ggsave(file.path("figures", "figure_8_sample_map.svg"),
       plot = my_plot, width = 10, height = 8, dpi = 900)
ggsave(file.path("figures", "figure_8_sample_map.pdf"),
       plot = my_plot, width = 10, height = 8, dpi = 900)
ggsave(file.path("figures", "figure_8_sample_map.png"),
       plot = my_plot, width = 10, height = 8, dpi = 900,
       bg = "white")

my_plot <- ggplot() +
  geom_sf(data = nyc_boroughs, fill = "lightgray", color = "black") +
  geom_sf(data = points_sf, aes(color = as.factor(F_bar)),
          size = 1.8, alpha = 1.0) +
  scale_color_manual(values = color_palette) +  # Customize the colors as needed
  theme_minimal() +
  labs(title = "Figure 8: Location of Sample Properties by FAR Limit",
       color = "FAR Limit")

ggsave(file.path("figures", "figure_8_sample_map_opaque.svg"),
       plot = my_plot, width = 10, height = 8, dpi = 900)
ggsave(file.path("figures", "figure_8_sample_map_opaque.pdf"),
       plot = my_plot, width = 10, height = 8, dpi = 900)
ggsave(file.path("figures", "figure_8_sample_map_opaque.png"),
       plot = my_plot, width = 10, height = 8, dpi = 900,
       bg = "white")