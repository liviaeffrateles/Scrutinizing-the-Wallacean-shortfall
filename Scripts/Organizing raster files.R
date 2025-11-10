library(terra)
library(dplyr)
library(raster)
library(mapview)
library(sf)

# This script organizes the rasters of the biogeographic realms
# that will be used in the environmental analyses

# Afrotropical ------------------------------------------------------------

# Import the 19 CHELSA bioclimatic variables aggregated to a 1° resolution
folder_path <- "D:/Scrutinizing-the-Wallacean-shortfall/Chelsa"
files <- list.files(folder_path, pattern = "^chelsaBio.*\\.tif$", full.names = TRUE)
chelsa_layers <- rast(files)

# Check the raster resolution
res(chelsa_layers)

# Crop the rasters to the Afrotropical realm
# Import the shapefile of biogeographic realms
BioRealms <- st_read("Grids/bioRealms.gpkg")

# Filter the Afrotropical realm
Afrotropical <- BioRealms %>%
  dplyr::filter(REALM == "AT")

# Convert to raster format
Afrotropical_vect <- terra::vect(Afrotropical)
Chelsa_Crop <- terra::crop(chelsa_layers, Afrotropical_vect)
climate_Afrotropical <- mask(Chelsa_Crop, Afrotropical_vect)
plot(climate_Afrotropical[["CHELSA_bio10_1981-2010_V.2.1"]])

# Get the coordinates of the bioclimatic variables
coords_AT <- as.data.frame(climate_Afrotropical, xy = T)

coords_AT <- coords_AT %>%
  mutate(y = round(y, 3),
         x = round(x, 3))

# Load Afrotropical grid
AT_grid <- st_read("Grids/AT_realm.gpkg")

mapview(AT_grid)

# Get the centroids of the grid
centroids_AT <- st_centroid(AT_grid)

centroids_AT <- centroids_AT %>% 
  dplyr::mutate(Longitude = sf::st_coordinates(.)[,1],
                Latitude = sf::st_coordinates(.)[,2]) 

centroids_AT$geometry <- NULL

centroids_AT <- centroids_AT %>%
  mutate(Latitude = round(Latitude, 1),
         Longitude = round(Longitude, 1))

# Identify coordinates common to the bioclimatic variables and the grid
coords_AT <- coords_AT %>%
  mutate(match = paste(y, x) %in% paste(centroids_AT$Latitude, centroids_AT$Longitude))

sum(coords_AT$match) # 2025

coords_AT_NA <- coords_AT %>%
  filter(match == "FALSE")

coords_AT_NA <- coords_AT_NA %>%
  filter(!is.na('CHELSA_bio1_1981-2010_V.2.1'))

coords_AT$bio <- coords_AT$'CHELSA_bio1_1981-2010_V.2.1' %in% coords_AT_NA$'CHELSA_bio1_1981-2010_V.2.1'
sum(coords_AT$bio)

coords_AT_re <- coords_AT %>%
  filter(bio == "TRUE") %>%
  dplyr::select(x, y)

coords_AT_re <- st_as_sf (coords_AT_re,
                          coords = c("x", "y"), crs = 4326)

mapview(coords_AT_re) +
  mapview(AT_grid)

# Convert points to SpatialPoints
points_AT <- as(coords_AT_re, "Spatial")

points_AT <- vect(points_AT) 

# Apply mask
climate_Afrotropical <- mask(climate_Afrotropical, points_AT, inverse = TRUE)
plot(climate_Afrotropical[[1]])

# Save the cropped and aggregated variables
for (i in 1:nlyr(climate_Afrotropical)) {
  layer_name <- paste0("Chelsa/climate_Afrotropical", i, ".tif")
  writeRaster(climate_Afrotropical[[i]], layer_name, overwrite = TRUE)
}

# Neotropical -------------------------------------------------------------

# Filter the Neotropical realm
Neotropical <- BioRealms %>%
  dplyr::filter(REALM == "NT")

# Convert to raster format
Neotropical_vect <- terra::vect(Neotropical)
Chelsa_Crop <- terra::crop(chelsa_layers, Neotropical_vect)
climate_Neotropical <- mask(Chelsa_Crop, Neotropical_vect)
plot(climate_Neotropical[["CHELSA_bio10_1981-2010_V.2.1"]])

NT_grid <- st_read("Grids/NT_realm.gpkg")

mapview(NT_grid)

# Get the coordinates of the bioclimatic variables
coords_NT <- as.data.frame(climate_Neotropical, xy = T)

coords_NT <- coords_NT %>%
  mutate(y = round(y, 3),
         x = round(x, 3))

# Get the centroids of the grid 
centroids_NT <- st_centroid(NT_grid)

centroids_NT <- centroids_NT %>% 
  dplyr::mutate(Longitude = sf::st_coordinates(.)[,1],
                Latitude = sf::st_coordinates(.)[,2]) 

centroids_NT$geometry <- NULL

centroids_NT <- centroids_NT %>%
  mutate(Latitude = round(Latitude, 1),
         Longitude = round(Longitude, 1))

# Identify coordinates common to the bioclimatic variables and the grid
coords_NT <- coords_NT %>%
  mutate(match = paste(y, x) %in% paste(centroids_NT$Latitude, centroids_NT$Longitude))

sum(coords_NT$match) # 2031

coords_NT_NA <- coords_NT %>%
  filter(match == "FALSE")

coords_NT_NA <- coords_NT_NA %>%
  filter(!is.na('CHELSA_bio18_1981-2010_V.2.1'))

coords_NT$bio <- coords_NT$'CHELSA_bio18_1981-2010_V.2.1' %in% coords_NT_NA$'CHELSA_bio18_1981-2010_V.2.1'
sum(coords_NT$bio)

coords_NT_re <- coords_NT %>%
  filter(bio == "TRUE") %>%
  dplyr::select(x, y)

coords_NT_re <- st_as_sf (coords_NT_re,
                          coords = c("x", "y"), crs = 4326)

mapview(coords_NT_re) +
  mapview(NT_grid)

# Convert points to SpatialPoints
points_NT <- as(coords_NT_re, "Spatial")
points_NT <- vect(points_NT) 

# Apply mask
climate_Neotropical <- mask(climate_Neotropical, points_NT, inverse = TRUE)
plot(climate_Neotropical[[1]])

# Save the cropped and aggregated variables
for (i in 1:nlyr(climate_Neotropical)) {
  layer_name <- paste0("Chelsa/climate_Neotropical", i, ".tif")
  writeRaster(climate_Neotropical[[i]], layer_name, overwrite = TRUE)
}


# Nearctic ----------------------------------------------------------------

NA_grid <- st_read("Grids/NA_realm.gpkg")
mapview(NA_grid)

# Filter the Nearctic realm
Nearctic <- BioRealms %>%
  dplyr::filter(REALM == "NA")

# Convert to raster format
Nearctic_vect <- terra::vect(Nearctic)
Chelsa_Crop <- terra::crop(chelsa_layers, Nearctic_vect)
climate_Nearctic <- mask(Chelsa_Crop, Nearctic_vect)
plot(climate_Nearctic[["CHELSA_bio10_1981-2010_V.2.1"]])

# Get the coordinates of the bioclimatic variables 
coords_NA <- as.data.frame(climate_Nearctic, xy = T)

coords_NA <- coords_NA %>%
  mutate(y = round(y, 3),
         x = round(x, 3))

# Get the centroids of the grid 
centroids_NA <- st_centroid(NA_grid)

centroids_NA <- centroids_NA %>% 
  dplyr::mutate(Longitude = sf::st_coordinates(.)[,1],
                Latitude = sf::st_coordinates(.)[,2]) 

centroids_NA$geometry <- NULL

centroids_NA <- centroids_NA %>%
  mutate(Latitude = round(Latitude, 1),
         Longitude = round(Longitude, 1))

# Identify coordinates common to the bioclimatic variables and the grid
coords_NA <- coords_NA %>%
  mutate(match = paste(y, x) %in% paste(centroids_NA$Latitude, centroids_NA$Longitude))

sum(coords_NA$match) #2230

coords_NA_NA <- coords_NA %>%
  filter(match == "FALSE")

coords_NA_NA <- coords_NA_NA %>%
  filter(!is.na('CHELSA_bio18_1981-2010_V.2.1'))

coords_NA$bio <- coords_NA$'CHELSA_bio18_1981-2010_V.2.1' %in% coords_NA_NA$'CHELSA_bio18_1981-2010_V.2.1'
sum(coords_NA$bio)

coords_NA_re <- coords_NA %>%
  filter(bio == "TRUE") %>%
  dplyr::select(x, y)

coords_NA_re <- st_as_sf (coords_NA_re,
                              coords = c("x", "y"), crs = 4326)

mapview(coords_NA_re) +
  mapview(NA_grid)

# Convert points to SpatialPoints
points_NA <- as(coords_NA_re, "Spatial")
points_NA <- vect(points_NA) 

# Apply mask
climate_Nearctic <- mask(climate_Nearctic, points_NA, inverse = TRUE)
plot(climate_Nearctic[[1]])

# Save the cropped and aggregated variables
for (i in 1:nlyr(climate_Nearctic)) {
  layer_name <- paste0("Chelsa/climate_Nearctic", i, ".tif")
  writeRaster(climate_Nearctic[[i]], layer_name, overwrite = TRUE)
}

# Indomalayan --------------------------------------------------------------
IM_grid <- st_read("grids/IM_realm.gpkg")
mapview(IM_grid)

# Filter the Indomalayan realm
Indomalayan <- BioRealms %>%
  dplyr::filter(REALM == "IM")

# Convert to raster format
Indomalayan_vect <- terra::vect(Indomalayan)
Chelsa_Crop <- terra::crop(chelsa_layers, Indomalayan_vect)
climate_Indomalayan <- mask(Chelsa_Crop, Indomalayan_vect)
plot(climate_Indomalayan[["CHELSA_bio10_1981-2010_V.2.1"]])

# Get the coordinates of the bioclimatic variables 
coords_IM <- as.data.frame(climate_Indomalayan, xy = T)

coords_IM <- coords_IM %>%
  mutate(y = round(y, 3),
         x = round(x, 3))

# Get the centroids of the grid 
centroids_IM <- st_centroid(IM_grid)

centroids_IM <- centroids_IM %>% 
  dplyr::mutate(Longitude = sf::st_coordinates(.)[,1],
                Latitude = sf::st_coordinates(.)[,2]) 

centroids_IM$geometry <- NULL

centroids_IM <- centroids_IM %>%
  mutate(Latitude = round(Latitude, 1),
         Longitude = round(Longitude, 1))

# Identify coordinates common to the bioclimatic variables and the grid
coords_IM <- coords_IM %>%
  mutate(match = paste(y, x) %in% paste(centroids_IM$Latitude, centroids_IM$Longitude))

sum(coords_IM$match) # 1085

coords_IM_NA <- coords_IM %>%
  filter(match == "FALSE")

coords_IM_NA <- coords_IM_NA %>%
  filter(!is.na('CHELSA_bio18_1981-2010_V.2.1'))

coords_IM$bio <- coords_IM$'CHELSA_bio18_1981-2010_V.2.1' %in% coords_IM_NA$'CHELSA_bio18_1981-2010_V.2.1'
sum(coords_IM$bio)

coords_IM_re <- coords_IM %>%
  filter(bio == "TRUE") %>%
  dplyr::select(x, y)

coords_IM_re <- st_as_sf (coords_IM_re,
                          coords = c("x", "y"), crs = 4326)

mapview(coords_IM_re) +
  mapview(IM_grid)

# Convert points to SpatialPoints
points_IM <- as(coords_IM_re, "Spatial")
points_IM <- vect(points_IM) 

# Apply mask
climate_Indomalayan <- mask(climate_Indomalayan, points_IM, inverse = TRUE)
plot(climate_Indomalayan[[1]])

# Save the cropped and aggregated variables
for (i in 1:nlyr(climate_Indomalayan)) {
  layer_name <- paste0("Chelsa/climate_Indomalayan", i, ".tif")
  writeRaster(climate_Indomalayan[[i]], layer_name, overwrite = TRUE)
}


# Eastern Palearctic ------------------------------------------------------

EP_grid <- st_read("Grids/EP_realm.gpkg")
mapview(EP_grid)

# Filter the Eastern Palearctic realm
E_Palearctic <- BioRealms %>%
  dplyr::filter(REALM == "EP")

# Convert to raster format
E_Palearctic_vect <- terra::vect(E_Palearctic)
Chelsa_Crop <- terra::crop(chelsa_layers, E_Palearctic_vect)
climate_E_Palearctic <- mask(Chelsa_Crop, E_Palearctic_vect)
plot(climate_E_Palearctic[["CHELSA_bio10_1981-2010_V.2.1"]])

# Get the coordinates of the bioclimatic variables 
coords_EP <- as.data.frame(climate_E_Palearctic, xy = T)

coords_EP <- coords_EP %>%
  mutate(y = round(y, 3),
         x = round(x, 3))

# Get the centroids of the grid 
centroids_EP <- st_centroid(EP_grid)

centroids_EP <- centroids_EP %>% 
  dplyr::mutate(Longitude = sf::st_coordinates(.)[,1],
                Latitude = sf::st_coordinates(.)[,2]) 

centroids_EP$geometry <- NULL

centroids_EP <- centroids_EP %>%
  mutate(Latitude = round(Latitude, 1),
         Longitude = round(Longitude, 1))

# Identify coordinates common to the bioclimatic variables and the grid
coords_EP <- coords_EP %>%
  mutate(match = paste(y, x) %in% paste(centroids_EP$Latitude, centroids_EP$Longitude))

sum(coords_EP$match) # 3612

coords_EP_NA <- coords_EP %>%
  filter(match == "FALSE")

coords_EP_NA <- coords_EP_NA %>%
  filter(!is.na('CHELSA_bio18_1981-2010_V.2.1'))

coords_EP$bio <- coords_EP$'CHELSA_bio18_1981-2010_V.2.1' %in% coords_EP_NA$'CHELSA_bio18_1981-2010_V.2.1'
sum(coords_EP$bio)

coords_EP_re <- coords_EP %>%
  filter(bio == "TRUE") %>%
  dplyr::select(x, y)

coords_EP_re <- st_as_sf (coords_EP_re,
                          coords = c("x", "y"), crs = 4326)

mapview(coords_EP_re) +
  mapview(EP_grid)

# Convert points to SpatialPoints
points_EP <- as(coords_EP_re, "Spatial")
points_EP <- vect(points_EP) 

# Apply mask
climate_E_Palearctic <- mask(climate_E_Palearctic, points_EP, inverse = TRUE)
plot(climate_E_Palearctic[[1]])

# Save the cropped and aggregated variables
for (i in 1:nlyr(climate_E_Palearctic)) {
  layer_name <- paste0("Chelsa/climate_E_Palearctic", i, ".tif")
  writeRaster(climate_E_Palearctic[[i]], layer_name, overwrite = TRUE)
}

# Western Palearctic ------------------------------------------------------

WP_grid <- st_read("grids/WP_realm.gpkg")
mapview(WP_grid)

# Filter the Western Palearctic realm
W_Palearctic <- BioRealms %>%
  dplyr::filter(REALM == "WP")

# Convert to raster format
W_Palearctic_vect <- terra::vect(W_Palearctic)
Chelsa_Crop <- terra::crop(chelsa_layers, W_Palearctic_vect)
climate_W_Palearctic <- mask(Chelsa_Crop, W_Palearctic_vect)
plot(climate_W_Palearctic[["CHELSA_bio10_1981-2010_V.2.1"]])

# Get the coordinates of the bioclimatic variables
coords_WP <- as.data.frame(climate_W_Palearctic, xy = T)

coords_WP <- coords_WP %>%
  mutate(y = round(y, 3),
         x = round(x, 3))

# Get the centroids of the grid
centroids_WP <- st_centroid(WP_grid)

centroids_WP <- centroids_WP %>% 
  dplyr::mutate(Longitude = sf::st_coordinates(.)[,1],
                Latitude = sf::st_coordinates(.)[,2]) 

centroids_WP$geometry <- NULL

centroids_WP <- centroids_WP %>%
  mutate(Latitude = round(Latitude, 1),
         Longitude = round(Longitude, 1))

# Identify coordinates common to the bioclimatic variables and the grid
coords_WP <- coords_WP %>%
  mutate(match = paste(y, x) %in% paste(centroids_WP$Latitude, centroids_WP$Longitude))

sum(coords_WP$match) # 2922

coords_WP_NA <- coords_WP %>%
  filter(match == "FALSE")

coords_WP_NA <- coords_WP_NA %>%
  filter(!is.na('CHELSA_bio18_1981-2010_V.2.1'))

coords_WP$bio <- coords_WP$'CHELSA_bio18_1981-2010_V.2.1' %in% coords_WP_NA$'CHELSA_bio18_1981-2010_V.2.1'
sum(coords_WP$bio)

coords_WP_re <- coords_WP %>%
  filter(bio == "TRUE") %>%
  dplyr::select(x, y)

coords_WP_re <- st_as_sf (coords_WP_re,
                          coords = c("x", "y"), crs = 4326)

mapview(coords_WP_re) +
  mapview(WP_grid)

# Convert points to SpatialPoints
points_WP <- as(coords_WP_re, "Spatial")
points_WP <- vect(points_WP) 

# Apply mask
climate_W_Palearctic <- mask(climate_W_Palearctic, points_WP, inverse = TRUE)
plot(climate_W_Palearctic[[1]])

# Save the cropped and aggregated variables
for (i in 1:nlyr(climate_W_Palearctic)) {
  layer_name <- paste0("Chelsa/climate_W_Palearctic", i, ".tif")
  writeRaster(climate_W_Palearctic[[i]], layer_name, overwrite = TRUE)
}


# Australasia -------------------------------------------------------------

AA_grid <- st_read("grids/AA_realm.gpkg")
mapview(AA_grid)

# Filter the Australasia realm
Australasia <- BioRealms %>%
  dplyr::filter(REALM == "AA")

# Convert to raster format
Australasia_vect <- terra::vect(Australasia)
Chelsa_Crop <- terra::crop(chelsa_layers, Australasia)
climate_Australasia <- mask(Chelsa_Crop, Australasia)
plot(climate_Australasia[["CHELSA_bio10_1981-2010_V.2.1"]])

# Get the coordinates of the bioclimatic variables
coords_AA <- as.data.frame(climate_Australasia, xy = T)

coords_AA <- coords_AA %>%
  mutate(y = round(y, 3),
         x = round(x, 3))


# Get the centroids of the grid
centroids_AA <- st_centroid(AA_grid)

centroids_AA <- centroids_AA %>% 
  dplyr::mutate(Longitude = sf::st_coordinates(.)[,1],
                Latitude = sf::st_coordinates(.)[,2]) 

centroids_AA$geometry <- NULL

centroids_AA <- centroids_AA %>%
  mutate(Latitude = round(Latitude, 1),
         Longitude = round(Longitude, 1))

# Identify coordinates common to the bioclimatic variables and the grid
coords_AA <- coords_AA %>%
  mutate(match = paste(y, x) %in% paste(centroids_AA$Latitude, centroids_AA$Longitude))

sum(coords_AA$match) # 1176

coords_AA_NA <- coords_AA %>%
  filter(match == "FALSE")

coords_AA_NA <- coords_AA_NA %>%
  filter(!is.na('CHELSA_bio18_1981-2010_V.2.1'))

coords_AA$bio <- coords_AA$'CHELSA_bio18_1981-2010_V.2.1' %in% coords_AA_NA$'CHELSA_bio18_1981-2010_V.2.1'
sum(coords_AA$bio)

coords_AA_re <- coords_AA %>%
  filter(bio == "TRUE") %>%
  dplyr::select(x, y)

coords_AA_re <- st_as_sf (coords_AA_re,
                          coords = c("x", "y"), crs = 4326)

mapview(coords_AA_re) +
  mapview(AA_grid)

# Convert points to SpatialPoints
points_AA <- as(coords_AA_re, "Spatial")
points_AA <- vect(points_AA) 

# Apply mask
climate_Australasia <- mask(climate_Australasia, points_AA, inverse = TRUE)
plot(climate_Australasia[[1]])

# Save the cropped and aggregated variables
for (i in 1:nlyr(climate_Australasia)) {
  layer_name <- paste0("Chelsa/climate_Australasia", i, ".tif")
  writeRaster(climate_Australasia[[i]], layer_name, overwrite = TRUE)
}
