library(sf)
library(dplyr)
library(rgdal)
library(tidyverse)
library(mapview)
library(KnowBR)
library(ggpubr)


# Inventory completeness analysis -----------------------------------------

# Load the data set
DataGBIF <- readRDS(here::here ("Data/RecordsGBIF.rds")) # GBIF
colnames(DataGBIF)
Type_locality <- readRDS("Data/Type_localities.rds") # Type_locality
colnames(Type_locality)

# Prepare the data set for analysis
DataGBIFKnow <- DataGBIF %>%
  dplyr::select(-scientificName, -Species_clean, -family, -year, -country, -country_suggested) %>%
  rename (Species = accepted_name) %>%
  rename(Longitude = decimalLongitude) %>%
  rename(Latitude = decimalLatitude) %>%
  dplyr::select(Species, Longitude, Latitude) 

colnames(DataGBIFKnow)

Type_localityKnow <- Type_locality %>%
  dplyr::select(-Type_locality, -Country, -REALM) %>%
  dplyr::select(Species, Longitude, Latitude) 

colnames(Type_localityKnow)

# Create a new data set with occurrence records from GBIF and type localities
DataGBIFKnowType <- bind_rows(DataGBIFKnow, Type_localityKnow)

DataGBIFKnowType <- DataGBIFKnowType %>%
  group_by(Species, Longitude, Latitude) %>%
  mutate(Count = n()) %>%
  ungroup() %>%
  distinct(Species, Longitude, Latitude, Count) %>%
  rename(Counts = Count)

DataGBIFKnowType <- DataGBIFKnowType %>%
  arrange(Species)

DataGBIFKnowType <- as.data.frame(DataGBIFKnowType)

# Convert the coordinates to an sf object 

DataGBIFKnowBR_sf <- st_as_sf(DataGBIFKnowType,
                              coords = c("Longitude", "Latitude"),
                              crs = 4326)

# Reproject the coordinates to an equal-area projection (Behrmann)
DataBehrmann_sf <- st_transform(DataGBIFKnowBR_sf, 
                                crs = "+proj=cea +lat_ts=30 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs") #World_Behrmann

# Extract longitude and latitude
DataBehrmann_sf <- DataBehrmann_sf %>% 
  dplyr::mutate(Longitude = sf::st_coordinates(.)[,1],
                Latitude = sf::st_coordinates(.)[,2])


DataBehrmann_sf$geometry <- NULL

BehrmannSnakesKnowBR_df <- DataBehrmann_sf %>%
  dplyr::select(Species, Longitude, Latitude, Counts)


# KnowBR -----------------------------------------------------------------

# Load user-defined grid
gridBehrmann_110km <- readOGR("Grids/global_grid.gpkg")

data(adworld)

# Calculate inventory completeness
RecordsGBIF_results <- KnowBPolygon(data = BehrmannSnakesKnowBR_df, 
                                    shape = gridBehrmann_110km, admAreas = FALSE,  # Use predefined grid as personalized polygons  
                                    shapenames = "id",  # WRITE HERE your unique "id" cell from the grid
                                    jpg = TRUE, dec = ".")

# Results -----------------------------------------------------------------

# Estimators from the inventory completeness analysis 
estimatorsBehrmann <- read.csv('Estimators.CSV', header = TRUE, sep = ",")

# Load user-defined grid
gridBehrmann_sf <- st_read("Grids/global_grid.gpkg")

# Rename the 'Area' column in the estimators to 'id' to match the 'id' column in the grid 
estimatorsBehrmann <- estimatorsBehrmann %>%
  rename (id = Area)

# Join inventory completeness estimators with the defined grid 
gridBehrmann_sf <- gridBehrmann_sf %>%
  left_join(estimatorsBehrmann, by = "id")

# Remove cells with no occurrence records
gridBehrmann_sf_2 <- gridBehrmann_sf %>% 
  filter(!is.na(Records))

# Plot completeness results establishing threshold of well-sampled
plot <- function(var, x, xtitle){
  ggplot(estimatorsBehrmann) +
    geom_point(aes(var, Completeness), pch = 19, size = 0.4) +
    geom_vline(xintercept = x, col = "#D55E00", lwd = 1, lty = 2) +
    geom_hline(yintercept = 70, col = 'grey', lwd = 1, lty = 2) + # Here completeness threshold = 70%
    theme_minimal() + 
    ylab('') +
    xlab(xtitle) +
    theme(strip.text.y = element_blank(),
          axis.text = element_text(size = 6),
          axis.title = element_text(size = 6))
}

# Distribution of completeness values by other estimators:
# Plot records value vs completeness (threshold set in 10):
p1 <- plot(estimatorsBehrmann$Records, 50, 'Records') + ylab('Completeness')

# Plot ratio value vs completeness (threshold set min 5):
p2 <- plot(estimatorsBehrmann$Ratio, 2, 'Ratio')

# Plot slope value vs completeness (threshold set in 0.1) to 
# Eastern Palearctic, Indo-Malay, Afrotropical, Neotropical:  
p3 <- plot(estimatorsBehrmann$Slope, 0.1, 'Slope')

# Plot slope value vs completeness (threshold set in 0.05) to
# Nearctic, Western Palearctic, Australasia:  
p4 <- plot(estimatorsBehrmann$Slope, 0.05, 'Slope')

figure.S2 <- ggarrange(p1, p2, p3, p4, ncol = 4, labels = c("(a)", "(b)", "(c)", "(d)"),
                        font.label = list(size = 7, face = "plain", color ="black"))
figure.S2

# Nearctic, Western Palearctic, Australasia
# Filter est dataset based on your chosen thresholds selecting well-sampled only:
estimatorsWS_1 <- gridBehrmann_sf %>% filter(Records > 50) %>% 
  filter(Ratio > 2) %>% 
  filter(Slope < 0.05)

# Transform landmass to the desired projection before creating the grid

estimatorsWS_WGS_1 <- st_transform(estimatorsWS_1, crs = 4326 )
max(estimatorsWS_WGS_1$Completeness)

# Extract centroids of well-sampled cells for the environmental Space analysis
WS_cent_1 <- st_centroid(estimatorsWS_WGS_1)

WS_cent_1 <- WS_cent_1 %>% 
  dplyr::mutate(Longitude = sf::st_coordinates(.)[,1],
                Latitude = sf::st_coordinates(.)[,2]) %>% 
  dplyr::select(c(id, Longitude, Latitude))


# Eastern Palearctic, Indomalayan, Afrotropical, Neotropical
# Filter the dataset based on your chosen thresholds selecting well-sampled cells only:
estimatorsWS_2 <- gridBehrmann_sf %>% filter(Records > 50) %>% 
  filter(Ratio > 2) %>% 
  filter(Slope < 0.1)

# Transform landmass to the desired projection before creating the grid
estimatorsWS_WGS_2 <- st_transform(estimatorsWS_2, crs = 4326 )

# Extract centroids of well-sampled for the environmental Space analysis
WS_cent_2 <- st_centroid(estimatorsWS_WGS_2)

WS_cent_2 <- WS_cent_2 %>% 
  dplyr::mutate(Longitude = sf::st_coordinates(.)[,1],
                Latitude = sf::st_coordinates(.)[,2]) %>% 
  dplyr::select(c(id, Longitude, Latitude))

#Separate well-sampled cells by biogeographic realms: NA, WP, AA, NT, AT, EP, IM

#Load grids
path <- "D:/Scrutinizing the Wallacean shortfall/Grids"
realms <- list.files(path = path, pattern = "_realm.*\\.gpkg$", full.names = TRUE)
all_realms <- lapply(realms, st_read)
names(all_realms) <- tools::file_path_sans_ext(basename(realms))

#Join the grids
AA_union <- st_union (all_realms[[1]])
WP_union <- st_union (all_realms[[7]])
NA_union <- st_union (all_realms[[5]])
AT_union <- st_union (all_realms[[2]])
EP_union <- st_union (all_realms[[3]])
IM_union <- st_union (all_realms[[4]])
NT_union <- st_union (all_realms[[6]])


well_sampled_points <- function(ws_sf, realm_union,realm_name) {
  result <- st_intersection(ws_sf, realm_union) %>%
    mutate(Longitude = st_coordinates(.)[,1],
           Latitude = st_coordinates(.)[,2]) %>%
    select(id, Longitude, Latitude) %>%
    mutate(Realm = realm_name)
  st_geometry(result) <- NULL
  return(result)}

well_sampled_AA_1 <- well_sampled_points(ws_cent_1_sf, AA_union, "AA")
well_sampled_NA_1 <- well_sampled_points(ws_cent_1_sf, NA_union, "NA")
well_sampled_WP_1 <- well_sampled_points(ws_cent_1_sf, WP_union, "WP")

well_sampled_AT_1 <- well_sampled_points(ws_cent_2_sf, AT_union, "AT")
well_sampled_NT_1 <- well_sampled_points(ws_cent_2_sf, NT_union, "NT")
well_sampled_IM_1 <- well_sampled_points(ws_cent_2_sf, IM_union, "IM")
well_sampled_EP_1 <- well_sampled_points(ws_cent_2_sf, EP_union, "EP")

well_sampled_1 <- bind_rows(well_sampled_AA_1, well_sampled_NA_1, well_sampled_WP_1)
well_sampled_2 <- bind_rows(well_sampled_AT_1, well_sampled_NT_1, well_sampled_IM_1, well_sampled_EP_1)

# well_sampled_1 <- write.csv("Processed_Data/well_sampled_1.csv") #save to use in other analyses 

# well_sampled_2 <- write.csv("Processed_Data/well_sampled_2.csv") #save to use in other analyses
