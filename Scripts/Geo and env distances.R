library(biscale)
library(sf)
library(openxlsx)
library(dplyr)
library(cowplot)
library(geosphere)
library(ggplot2)
library(rnaturalearth)


# This script contains the code to calculate geographic and environmental distances, 
# generate their respective figures, and create bivariate maps. 


# Geographic distances -----------------------------------------------------

calculate_geographic_distances <- function(realm_code, grid_path, well_sampled_data) {
  grid <- st_read(grid_path, quiet = TRUE)
  grid_proj <- st_transform(grid, crs = "+proj=cea +lat_ts=30 +lon_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")
  
  centroids <- st_centroid(grid_proj)
  centroids_WGS <- st_transform(centroids, crs = 4326)
  
  centroids_df <- centroids_WGS %>%
    mutate(Longitude = st_coordinates(.)[, 1],
           Latitude = st_coordinates(.)[, 2],
           id = 1:n()) %>%
    select(id, Longitude, Latitude) %>%
    st_drop_geometry()
  
  WS_filtered <- well_sampled_data %>%
    filter(Realm == realm_code) %>%
    mutate(id_WS = row_number()) %>%
    select(id_WS, Longitude_WS = Longitude, Latitude_WS = Latitude)
  
  distances <- expand.grid(id = centroids_df$id, id_WS = WS_filtered$id_WS) %>%
    left_join(centroids_df, by = "id") %>%
    left_join(WS_filtered, by = "id_WS") %>%
    rowwise() %>%
    mutate(distance_km = if (any(is.na(c(Longitude, Latitude, Longitude_WS, Latitude_WS)))) NA_real_
           else distGeo(c(Longitude, Latitude), c(Longitude_WS, Latitude_WS)) / 1000) %>%
    ungroup()
  
  result <- distances %>%
    group_by(id) %>%
    slice_min(distance_km, n = 1, na_rm = TRUE) %>%
    ungroup() %>%
    mutate(distance_km = round(distance_km, 4),
           scaled_data = (distance_km - min(distance_km, na.rm = TRUE)) / 
             (max(distance_km, na.rm = TRUE) - min(distance_km, na.rm = TRUE)),
           realm = realm_code) %>%
    select(realm, id, Longitude, Latitude, distance_km, scaled_data)
  
  return(result)
}

# Load well-sampled data obtained from the inventory completeness analysis
# well_sampled_1 <- read.csv("Processed_Data/well_sampled_1.csv")
# well_sampled_2 <- read.csv("Processed_Data/well_sampled_2.csv")
well_sampled_1 <- well_sampled_1 %>%
  mutate(Realm = ifelse(is.na(Realm), "NA", Realm))

WS <- bind_rows(well_sampled_1, well_sampled_2)

# List of realm codes and shapefile paths
realms <- c("AA", "AT", "EP", "IM", "NA", "NT", "WP")
shapefile_paths <- paste0("Grids/", realms, "_realm.gpkg")
names(shapefile_paths) <- realms

# Apply the function to all realms
results <- lapply(realms, function(realm) {
  calculate_geographic_distances(realm, shapefile_paths[[realm]], WS)
})

all_geo_distances <- bind_rows(results)

#Filter geographic distances by realm

geo_dist_NT <- all_geo_distances %>% 
  filter(realm == "NT")
geo_dist_AT <- all_geo_distances %>%
  filter(realm == "AT")
geo_dist_IM <- all_geo_distances %>% 
  filter(realm == "IM")
geo_dist_IM <- geo_dist_IM %>% distinct()
geo_dist_EP <- all_geo_distances %>% 
  filter(realm == "EP")
geo_dist_WP <- all_geo_distances %>% 
  filter(realm == "WP")
geo_dist_AA <- all_geo_distances %>% 
  filter(realm == "AA")
geo_dist_NA <- all_geo_distances %>% 
  filter(realm == "NA")
geo_dist_NA <- geo_dist_NA %>% distinct()

# Filter well-sampled cells for each realm

WS_NT <- WS %>% 
  filter(Realm == "NT")
WS_AT <- WS %>% 
  filter(Realm == "AT")
WS_IM <- WS %>% 
  filter(Realm == "IM")
WS_EP <- WS %>% 
  filter(Realm == "EP")
WS_WP <- WS %>% 
  filter(Realm == "WP")
WS_AA <- WS %>% 
  filter(Realm == "AA")
WS_NA <- WS %>% 
  filter(Realm == "NA") 

# Figure - Geo distance --------------------------------------------------------

world_sf <- ne_countries(scale = "medium", returnclass = "sf") %>%
  dplyr::filter(name != "Antarctica")

cols_geo <- colorRampPalette(c("#88CCEE", "#44AA99", "#DDCC77", "#CC6677"))

geo_distances <- function(geo_distance, ws_data, xlim, ylim, legend_pos, legend_dir, guide_width, guide_height, 
                          legend_hjust) {
  ggplot() + 
    geom_tile(data = geo_distance, aes(x = Longitude, y = Latitude, fill = distance_km)) +
    geom_sf(data = world_sf, fill = NA, color = "black") +
    geom_tile(data = ws_data, aes(x = Longitude, y = Latitude), width = 1, height = 1, linewidth = 0.2, color = "black", fill = NA) +
    coord_sf(ylim = ylim, xlim = xlim) +
    theme_test(base_size = 10) +
    labs(fill = "Geographic
distance", x = "Longitude", y = "Latitude") +
    theme(panel.background = element_rect(fill = "#F7FBFF"),
          legend.background = element_rect(fill = NA, color = NA), 
          legend.key = element_rect(fill = NA, color = NA)) +
    ggplot2::theme(
      plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "mm"),
      legend.text = element_text(size = 8), 
      legend.title = element_text(size = 8),
      axis.text = element_text(size = 10),
      axis.title.x = element_text(size = 10),
      axis.title.y = element_text(size = 10),
      axis.ticks.y = element_line(),
      legend.key.height = unit(0.4, 'cm'),
      legend.key.width = unit(0.4, 'cm'),
      legend.position = "inside",
      legend.position.inside = legend_pos, 
      legend.direction = legend_dir) +
    guides(fill = guide_colorbar(title.position = "top", title.hjust = legend_hjust, 
                                 barwidth = guide_width, barheight = guide_height)) +
    scale_fill_gradientn(colors = cols_geo(10))}

geo_IM_map <- geo_distances(geo_dist_IM, WS_IM, c(64, 135), c(-12.5, 36),
                            legend_pos = c(0.3, .13), legend_dir = "horizontal",
                            guide_width = 5, guide_height = 0.5, legend_hjust = 0.5)

geo_AT_map <- geo_distances(geo_dist_AT, WS_AT, c(-18, 78), c(-40.5, 25.9),
                            legend_pos = c(0.15, 0.3), legend_dir = "vertical",
                            guide_width = 0.5, guide_height = 4, legend_hjust = 0)

geo_NT_map <- geo_distances(geo_dist_NT, WS_NT, c(-140, -7), c(-58.5, 31.9),
                            legend_pos = c(0.1, .25), legend_dir = "vertical",
                            guide_width = 0.5, guide_height = 4, legend_hjust = 0)

geo_EP_map <- geo_distances(geo_dist_EP, WS_EP, c(45, 174), c(22, 79),
                            legend_pos = c(0.8, .13), legend_dir = "vertical",
                            guide_width = 0.5, guide_height = 4, legend_hjust = 0)

geo_NA_map <- geo_distances(geo_dist_NA, WS_NA, c(-167.2, -18), c(20, 83.6),
                            legend_pos = c(0.1, .25), legend_dir = "vertical",
                            guide_width = 0.5, guide_height = 4, legend_hjust = 0)

geo_AA_map <- geo_distances(geo_dist_AA, WS_AA, c(75.7, 180), c(-58, 7.5),
                            legend_pos = c(0.1, .25), legend_dir = "vertical",
                            guide_width = 0.5, guide_height = 4, legend_hjust = 0)

geo_WP_map <- geo_distances(geo_dist_WP, WS_WP, c(-47, 90), c(16, 80),
                            legend_pos = c(0.1, .25), legend_dir = "vertical",
                            guide_width = 0.5, guide_height = 4, legend_hjust = 0)


# Environmental distances --------------------------------------------------

calculate_environmental_distance <- function(realm_code) {
  # Load input files obtained from environmental analysis 
  # ID_WS <- read.table(paste0("Processed_Data/values_WS_", realm_code, ".txt"))
  # coords_WS <- read.table(paste0("Processed_Data/coords_WS_", realm_code, ".txt"))
  # coords_PCA <- read.table(paste0("Processed_Data/coords_PCA_", realm_code, ".txt"))
  
  # Prepare PCA_WS: bind IDs and PCA, join coordinates, assign id_WS
  PCA_WS <- bind_cols(as.data.frame(ID_WS), as.data.frame(coords_WS)) %>%
    rename(id = x) %>%
    left_join(coords_PCA %>% mutate(id = 1:n()), by = "id") %>%
    mutate(id_WS = 1:n())
  
  PCA_WS_2 <- PCA_WS %>%
    select(id_WS, PC1.x, PC2.x) %>%
    rename(PC1_WS = PC1.x, PC2_WS = PC2.x)
  
  # Prepare coords_PCA_2: drop NA and reformat
  coords_PCA_2 <- coords_PCA %>%
    filter(!is.na(PC1)) %>%
    select(PC1, PC2) %>%
    mutate(id_AT = 1:n()) %>%
    rename(PC1_AT = PC1, PC2_AT = PC2)
  
  # Calculate environmental distances (Euclidean)
  distances_env <- expand.grid(id_AT = coords_PCA_2$id_AT, id_WS = PCA_WS_2$id_WS) %>%
    left_join(coords_PCA_2, by = "id_AT") %>%
    left_join(PCA_WS_2, by = "id_WS") %>%
    rowwise() %>%
    mutate(distance_env = sqrt((PC1_AT - PC1_WS)^2 + (PC2_AT - PC2_WS)^2)) %>%
    ungroup()
  
  # Keep shortest distance per grid cell
  distances_env_s <- distances_env %>%
    group_by(id_AT) %>%
    slice_min(distance_env, n = 1) %>%
    ungroup()
  
  # Join final coordinates
  coords_PCA_3 <- coords_PCA %>%
    filter(!is.na(PC1)) %>%
    select(x, y, PC1, PC2)
  
  result <- distances_env_s %>%
    left_join(coords_PCA_3, by = c("PC1_AT" = "PC1")) %>%
    mutate(realm = realm_code)
  
  return(result)
}

realms <- c("NT", "AT", "IM", "EP", "WP", "AA", "NA")

env_distances <- lapply(realms, calculate_environmental_distance)

all_env_distances <- bind_rows(env_distances)

# Filter environmental distances for each realm

env_dist_NT <- all_env_distances %>% 
  filter(realm == "NT")
env_dist_AT <- all_env_distances %>%
  filter(realm == "AT")
env_dist_IM <- all_env_distances %>% 
  filter(realm == "IM")
env_dist_EP <- all_env_distances %>% 
  filter(realm == "EP")
env_dist_WP <- all_env_distances %>% 
  filter(realm == "WP")
env_dist_AA <- all_env_distances %>% 
  filter(realm == "AA")
env_dist_NA <- all_env_distances %>% 
  filter(realm == "NA")

# Figure - Environmental distances ----------------------------------------

cols_env <- colorRampPalette(c("#56B4E9", "#70C079", "#E69F00", "#CC79A7"))

env_distances <- function (env_distance, ws_data, xlim, ylim, legend_pos, legend_dir, guide_width, guide_height, 
                           legend_hjust) {
  ggplot() +
    geom_tile(data = env_distance, aes(x = x, y = y, fill = distance_env)) +
    geom_tile(data = ws_data, aes(x = Longitude, y = Latitude), width = 1, height = 1, linewidth = 0.2, color = "black", fill = NA) +
    scale_fill_gradientn(colors = cols_env(10)) +
    geom_sf(data = world_sf, fill = NA, color = "black") +
    coord_sf(ylim = ylim, xlim = xlim) +
    theme_test(base_size = 10) +
    theme(panel.background = element_rect(fill = "#F7FBFF"),
          legend.background = element_rect(fill = NA, color = NA), 
          legend.key = element_rect(fill = NA, color = NA))  +
    ggplot2::theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "mm"),
                   legend.text = element_text(size = 8), 
                   legend.title = element_text(size = 8),
                   axis.text = element_text(size = 10),
                   axis.title.x = element_text(size = 10),
                   axis.title.y = element_text(size = 10),
                   axis.ticks.y = element_line(),
                   legend.key.height = unit(0.4, 'cm'),
                   legend.key.width = unit(0.4, 'cm'),
                   legend.position = "inside",
                   legend.position.inside = legend_pos,
                   legend.direction = legend_dir) +
    guides(fill = guide_colorbar(title.position = "top", title.hjust = legend_hjust, 
                                 barwidth = guide_width, barheight = guide_height)) +
    labs(fill = "Environmental 
distance", x = "Longitude", y = "Latitude")}

env_IM_map <- env_distances(env_dist_IM, WS_IM, c(64, 135), c(-12.5, 36),
                            legend_pos = c(0.3, .13), legend_dir = "horizontal",
                            guide_width = 4, guide_height = 0.5, legend_hjust = 0.5)

env_AT_map <- env_distances(env_dist_AT, WS_AT, c(-18, 78), c(-40.5, 25.9),
                            legend_pos = c(0.15, 0.3), legend_dir = "vertical",
                            guide_width = 0.5, guide_height = 4, legend_hjust = 0)

env_NT_map <- env_distances(env_dist_NT, WS_NT, c(-140, -7), c(-58.5, 31.9),
                            legend_pos = c(0.15, .28), legend_dir = "vertical",
                            guide_width = 0.5, guide_height = 4, legend_hjust = 0)

env_EP_map <- env_distances(env_dist_EP, WS_EP, c(45, 174), c(22, 79),
                            legend_pos = c(0.83, .13), legend_dir = "horizontal",
                            guide_width = 5, guide_height = 0.5, legend_hjust = 0.5)

env_NA_map <- env_distances(env_dist_NA, WS_NA, c(-167.2, -18), c(20, 83.6),
                            legend_pos = c(0.15, .28), legend_dir = "vertical",
                            guide_width = 0.5, guide_height = 4, legend_hjust = 0)

env_AA_map <- env_distances(env_dist_AA, WS_AA, c(75.7, 180), c(-58, 7.5),
                            legend_pos = c(0.15, .28), legend_dir = "vertical",
                            guide_width = 0.5, guide_height = 4, legend_hjust = 0)

env_WP_map <- env_distances(env_dist_WP, WS_WP, c(-47, 90), c(16, 80),
                            legend_pos = c(0.15, .28), legend_dir = "vertical",
                            guide_width = 0.5, guide_height = 4, legend_hjust = 0)



# Bivariate maps ----------------------------------------------------------

# Join the previously calculated geographic and environmental distances to the 
# grids of each biogeographic realm

process_grid <- function(realm_code, geo_dist, env_dist, dec_lat = 2, dec_long = 2) {
  
  # Load grid shapefiles for each biogeographic realm
  grid_path <- paste0("Grids/", realm_code, "_realm.gpkg")
  grid <- st_read(grid_path, quiet = TRUE)
  
  # Adjust geographic distances by setting values < 0.5 km to 0
  geo_dist <- geo_dist %>%
    mutate(dist_km = ifelse(distance_km < 0.5, 0, distance_km))
  
  # Join geographic distances with the grid
  grid_cent <- bind_cols(grid, geo_dist)
  
  # Round geographic coordinates
  grid_cent <- grid_cent %>%
    mutate(Longitude = round(Longitude, dec_long),
           Latitude = round(Latitude, dec_lat))
  
  # Adjust environmental distances by setting values < 0.001 to 0
  env_dist <- env_dist %>%
    mutate(dist_env = ifelse(distance_env < 0.001, 0, distance_env),
           x = round(x, dec_long),
           y = round(y, dec_lat))
  
  # Join data frames by matching geographic coordinates
  grid_cent <- grid_cent %>%
    inner_join(env_dist, by = c("Longitude" = "x", "Latitude" = "y"))
  
  # Select specific columns
  grid_cent %>%
    select(distance_km, dist_km, distance_env, dist_env, realm.x, geom) %>%
    rename(realm = realm.x)
}

# Apply the function to each realm

NT_grid_cent <- process_grid("NT", geo_dist_NT, env_dist_NT)
IM_grid_cent <- process_grid("IM", geo_dist_IM, env_dist_IM)
AT_grid_cent <- process_grid("AT", geo_dist_AT, env_dist_AT)
EP_grid_cent <- process_grid("EP", geo_dist_EP, env_dist_EP)
WP_grid_cent <- process_grid("WP", geo_dist_WP, env_dist_WP, dec_lat = 1)
NA_grid_cent <- process_grid("NA", geo_dist_NA, env_dist_NA)
AA_grid_cent <- process_grid("AA", geo_dist_AA, env_dist_AA)

# Create classes for each biogeographic realm

# NT
data_NT <- bi_class(NT_grid_cent, x = dist_km, y = dist_env, 
                    style = "quantile", dim = 4)

# AA
data_AA <- bi_class(AA_grid_cent, x = dist_km, y = dist_env, 
                    style = "quantile", dim = 4)

#WP
data_WP <- bi_class(WP_grid_cent, x = dist_km, y = dist_env, 
                    style = "quantile", dim = 4)

# EP
data_EP <- bi_class(EP_grid_cent, x = dist_km, y = dist_env, 
                    style = "quantile", dim = 4)

# AT
data_AT <- bi_class(AT_grid_cent, x = dist_km, y = dist_env, 
                    style = "quantile", dim = 4)

# IM 
data_IM <- bi_class(IM_grid_cent, x = dist_km, y = dist_env, 
                    style = "quantile", dim = 4)

# NA
NA_grid_cent <- NA_grid_cent %>%
  mutate(dist_km = dist_km + runif(n(), -0.0000001, 0.0000001),
         dist_env = dist_env + runif(n(), -0.0000001, 0.0000001))

data_NA <- bi_class(NA_grid_cent, x = dist_km, y = dist_env, 
                    style = "quantile", dim = 4)

data_NA <- data_NA %>%
  mutate(bi_class = case_when(
    distance_km < 0.5 ~ "1-1",
    distance_env < 0.001 ~ "1-1",
    TRUE ~ bi_class
  ))



# Figure - bivariate maps -------------------------------------------------


bivariate_map <- function(bi_data, ws_data, xlim, ylim, lwidth) {
  ggplot() +
  geom_sf(data = bi_data, color = NA, mapping = aes(fill = bi_class), size = 0.1, show.legend = F) +
  geom_sf(data = world_sf, fill = NA) +
  geom_tile(data = ws_data, aes(x = Longitude, y = Latitude), width = 1, height = 1, linewidth = lwidth, color = "black", fill = NA) +
  coord_sf(ylim = ylim, xlim = xlim) +
  labs(x = "Longitude", y = "Latitude")+
  bi_scale_fill(pal = "BlueYl", dim = 4) +
  labs(
  ) +
  bi_theme()+
  theme_test(base_size = 10) +
  theme(panel.background = element_rect(fill = "#F7FBFF"),
        legend.background = element_rect(fill = NA, color = NA), 
        legend.key = element_rect(fill = NA, color = NA)) + 
  ggplot2::theme(
    plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "mm"),
    legend.text = element_text(size = 6.5), 
    axis.text = element_text(size = 10),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    legend.key.height = unit(0.4, 'cm'),
    legend.key.width = unit(0.4, 'cm'))}

legend <- bi_legend(pal = "BlueYl",
                     dim = 4,
                    xlab = "Geo dist",
                    ylab = "Env dist",
                    size = 6.5) +
  theme(panel.background = element_rect(fill = NA, color = NA), 
        plot.background = element_rect(fill = NA, color = NA))

bivariate_IM_map <- bivariate_map(data_IM, WS_IM, c(64, 135), c(-12.5, 36), lwidth = 0.2)

bivariate_AT_map <- bivariate_map(data_AT, WS_AT, c(-18, 78), c(-40.5, 25.9), lwidth = 0.2)

bivariate_NT_map <- bivariate_map(data_NT, WS_NT, c(-140, -7), c(-58.5, 31.9), lwidth = 0.2)

bivariate_EP_map <- bivariate_map(data_EP, WS_EP, c(45, 174), c(22, 79), lwidth = 0.2)

bivariate_NA_map <- bivariate_map(data_NA, WS_NA, c(-167.2, -18), c(20, 83.6), lwidth = 0.2)

bivariate_AA_map <- bivariate_map(data_AA, WS_AA, c(75.7, 180), c(-58, 7.5), lwidth = 0.2)

bivariate_WP_map <- bivariate_map(data_WP, WS_WP, c(-47, 90), c(16, 80), lwidth = 0.2)


finalPlot_IM <- ggdraw() +
  draw_plot(bivariate_IM_map, 0, 0, 1, 1) +
  draw_plot(legend, 0.2, .12, 0.29, 0.29)


finalPlot_AT <- ggdraw() +
  draw_plot(bivariate_AT_map, 0, 0, 1, 1) +
  draw_plot(legend, 0.15, .20, 0.29, 0.29)


finalPlot_NT <- ggdraw() +
  draw_plot(bivariate_NT_map, 0, 0, 1, 1) +
  draw_plot(legend, 0.18, .25, 0.29, 0.29)

finalPlot_EP <- ggdraw() +
  draw_plot(bivariate_EP_map, 0, 0, 1, 1) +
  draw_plot(legend, 0.6, .12, 0.55, 0.29)

finalPlot_NA <- ggdraw() +
  draw_plot(bivariate_NA_map, 0, 0, 1, 1) +
  draw_plot(legend, 0.15, .25, 0.29, 0.29)

finalPlot_AA <- ggdraw() +
  draw_plot(bivariate_AA_map, 0, 0, 1, 1) +
  draw_plot(legend, 0.15, .15, 0.29, 0.29)

finalPlot_WP <- ggdraw() +
  draw_plot(bivariate_WP_map, 0, 0, 1, 1) + 
  draw_plot(legend, 0.08, .20, 0.29, 0.29)

# Biogeographic realms

data_bivariate <- bind_rows(data_AA, data_AT, data_EP, data_IM, data_NA, data_NT, data_WP)

bivariate_global_map <- bivariate_map(data_bivariate, WS, c(-163.5, 163.5), c(-58, 75), lwidth = 0.12)

finalPlot_realms <- ggdraw() +
  draw_plot(bivariate_global_map, 0, 0, 1, 1) +
  draw_plot(legend, 0.05, .25, 0.29, 0.29)

