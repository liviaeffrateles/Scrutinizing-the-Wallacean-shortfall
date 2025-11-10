library(sf)
library(mapview)
library(dplyr)
library(ggplot2)
library(cowplot)
library(ggpubr)
library(patchwork)

# This script contains the code required to produce Figures S1 from the Supporting Information 
# and Figures 1 and 2 from the main text

# Biogeographic realms - Figure S1 ---------------------------------------------

bioRealms <- st_read("Grids/bioRealms.gpkg")

world_sf <- ne_countries(scale = "medium", returnclass = "sf") %>%
  dplyr::filter(name != "Antarctica")

cols <- colorRampPalette (c("#FB9A99", "#B2DF8A", "#A6CEE3", "#FDBF6F", "#1F78B4", "#FF7F00", "#E31A1C", "#DEEBF7")) 

realms_map <- ggplot() +
  geom_sf (data = world_sf, fill = NA, color = "black") +
  geom_sf(data = bioRealms, aes(fill = REALM)) +
  coord_sf(ylim = c(-58, 75), xlim = c(-163.5, 163.5)) +
  scale_fill_manual(values = cols(8), breaks = setdiff(unique(bioRealms$REALM), "OC")) +
  theme_test(base_size = 10) +
  theme(panel.background = element_rect(fill = "#F7FBFF"),
        legend.background = element_rect(fill = "#F7FBFF"), 
        legend.key = element_rect(fill = "#F7FBFF")) +
  ggplot2::theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "mm"),
                 legend.text = element_text(size = 6.5),
                 legend.title = element_text(size = 6.5), 
                 legend.key.height = unit(0.4, 'cm'),
                 legend.key.width = unit(0.4, 'cm'),
                 legend.position = c(0.1, .3),
                 legend.key.size = unit(0.5, "lines"),
                 axis.ticks = element_blank(),  
                 axis.text.x = element_blank(),  
                 axis.text.y = element_blank(),  
                 axis.title.x = element_blank(),  
                 axis.title.y = element_blank()) +
  labs(fill = "Biogeographic realm", x = "Longitude", y = "Latitude")

realms_map

# Histograms Figure 2 ----------------------------------------------------------

# Load the data set
DataGBIF <- readRDS(here::here ("Data/RecordsGBIF.rds")) #GBIF
colnames(DataGBIF)

# Load biogeographic realms
Realm <- st_read("Grids/bioRealms.gpkg")
mapview(Realm)

# Separate each biogeographic realm
NA_realm <- Realm %>%
  filter(REALM == "NA")

mapview(NA_realm)

NT_realm <- Realm %>%
  filter(REALM == "NT")

mapview(NT_realm)

WP_realm <- Realm %>%
  filter(REALM == "WP")

mapview(WP_realm)

EP_realm <- Realm %>%
  filter(REALM == "EP")

mapview(EP_realm)

AT_realm <- Realm %>%
  filter(REALM == "AT")

mapview(AT_realm)

IM_realm <- Realm %>%
  filter(REALM == "IM")

mapview(IM_realm)

AA_realm <- Realm %>%
  filter(REALM == "AA")

mapview(AA_realm)

OC_realm <- Realm %>%
  filter(REALM == "OC")

mapview(OC_realm)


# Convert coordinates to an sf object
DataGBIF_sf <- st_as_sf(DataGBIF,
                        coords = c("decimalLongitude", "decimalLatitude"),
                        crs = 4326)


# Group occurrence records by biogeographic realm 

#NA
Records_NA <- st_intersection(DataGBIF_sf, NA_realm) #614845

Records_NA <- Records_NA %>% 
  dplyr::mutate(Longitude = sf::st_coordinates(.)[,1],
                Latitude = sf::st_coordinates(.)[,2]) %>%
  select(accepted_name, family, year, REALM, Longitude,  Latitude, geometry) %>%
  arrange(year)

mapview(Records_NA) +
  mapview(NA_realm)

Records_NA$geometry <- NULL

# Count the number of occurrence records per species
NA_counts <- Records_NA %>%
  group_by(accepted_name) %>%
  summarise(n = n())

hist(NA_counts$n)

# NT
Records_NT <- st_intersection(DataGBIF_sf, NT_realm) #145716

Records_NT <- Records_NT %>% 
  dplyr::mutate(Longitude = sf::st_coordinates(.)[,1],
                Latitude = sf::st_coordinates(.)[,2]) %>%
  select(accepted_name, family, year, REALM, Longitude,  Latitude, geometry) %>%
  arrange(year)

mapview(Records_NT) +
  mapview(NT_realm)

Records_NT$geometry <- NULL


# Count the number of occurrence records per species
NT_counts <- Records_NT %>%
  group_by(accepted_name) %>%
  summarise(n = n())

hist(NT_counts$n)

# WP
Records_WP <- st_intersection(DataGBIF_sf, WP_realm) #361127

Records_WP <- Records_WP %>% 
  dplyr::mutate(Longitude = sf::st_coordinates(.)[,1],
                Latitude = sf::st_coordinates(.)[,2]) %>%
  select(accepted_name, family, year, REALM, Longitude,  Latitude, geometry) %>%
  arrange(year)

mapview(Records_WP) +
  mapview(WP_realm)

Records_WP$geometry <- NULL

# Count the number of occurrence records per species
WP_counts <- Records_WP %>%
  group_by(accepted_name) %>%
  summarise(n = n())

hist(WP_counts$n)

# EP
Records_EP <- st_intersection(DataGBIF_sf, EP_realm) #16383

Records_EP <- Records_EP %>% 
  dplyr::mutate(Longitude = sf::st_coordinates(.)[,1],
                Latitude = sf::st_coordinates(.)[,2]) %>%
  select(accepted_name, family, year, REALM, Longitude,  Latitude, geometry) %>%
  arrange(year)

mapview(Records_EP) +
  mapview(EP_realm)

Records_EP$geometry <- NULL

# Count the number of occurrence records per species
EP_counts <- Records_EP %>%
  group_by(accepted_name) %>%
  summarise(n = n())

hist(EP_counts$n)

# AT
Records_AT <- st_intersection(DataGBIF_sf, AT_realm) #42411

Records_AT <- Records_AT %>% 
  dplyr::mutate(Longitude = sf::st_coordinates(.)[,1],
                Latitude = sf::st_coordinates(.)[,2]) %>%
  select(accepted_name, family, year, REALM, Longitude,  Latitude, geometry) %>%
  arrange(year)

mapview(Records_AT) +
  mapview(AT_realm)

Records_AT$geometry <- NULL

# Count the number of occurrence records per species
AT_counts <- Records_AT %>%
  group_by(accepted_name) %>%
  summarise(n = n())

hist(AT_counts$n)

# IM
Records_IM <- st_intersection(DataGBIF_sf, IM_realm) #52479

Records_IM <- Records_IM %>% 
  dplyr::mutate(Longitude = sf::st_coordinates(.)[,1],
                Latitude = sf::st_coordinates(.)[,2]) %>%
  select(accepted_name, family, year, REALM, Longitude,  Latitude, geometry) %>%
  arrange(year)

mapview(Records_IM) +
  mapview(IM_realm)

Records_IM$geometry <- NULL

# Count the number of occurrence records per species
IM_counts <- Records_IM %>%
  group_by(accepted_name) %>%
  summarise(n = n())

hist(IM_counts$n)

# AA
Records_AA <- st_intersection(DataGBIF_sf, AA_realm) #133700

Records_AA <- Records_AA %>% 
  dplyr::mutate(Longitude = sf::st_coordinates(.)[,1],
                Latitude = sf::st_coordinates(.)[,2]) %>%
  select(accepted_name, family, year, REALM, Longitude,  Latitude, geometry) %>%
  arrange(year)

mapview(Records_AA) +
  mapview(AA_realm)

Records_AA$geometry <- NULL

# Count the number of occurrence records per species
AA_counts <- Records_AA %>%
  group_by(accepted_name) %>%
  summarise(n = n())

hist(AA_counts$n)

# OC
Records_OC <- st_intersection(DataGBIF_sf, OC_realm) #18

Records_OC <- Records_OC %>% 
  dplyr::mutate(Longitude = sf::st_coordinates(.)[,1],
                Latitude = sf::st_coordinates(.)[,2]) %>%
  select(accepted_name, family, year, REALM, Longitude,  Latitude, geometry) %>%
  arrange(year)

mapview(Records_OC) +
  mapview(OC_realm)

Records_OC$geometry <- NULL

# Count the number of occurrence records per species
OC_counts <- Records_OC %>%
  group_by(accepted_name) %>%
  summarise(n = n())

hist(OC_counts$n)


# Count the number of all occurrence records per species 
GBIF_counts <- DataGBIF %>%
  group_by(accepted_name) %>%
  summarise(n = n())

hist(GBIF_counts$n)

# Realms ------------------------------------------------------------------

realms <- st_read("Grids/realms.gpkg")
mapview(realms)

bio_NA <- realms %>%
  filter(REALM == "NA")

bio_NT <- realms %>%
  filter(REALM == "NT")

bio_AA <- realms %>%
  filter(REALM == "AA")

bio_EP <- realms %>%
  filter(REALM == "EP")

bio_WP <- realms %>%
  filter(REALM == "WP")

bio_IM <- realms %>%
  filter(REALM == "IM")

bio_AT <- realms %>%
  filter(REALM == "AT")

bio_NA_map <- ggplot() + 
  geom_sf(data = bio_NA, fill = "lightgrey", color = "lightgrey") +
  coord_sf(ylim = c(20, 85), xlim = c(-172, -18)) + 
  theme_void() + 
  theme(plot.background = element_rect(fill = "transparent", color = NA))

bio_NA_map

bio_NT_map <- ggplot() + 
  geom_sf(data = bio_NT, fill = "lightgrey", color = "lightgrey") +
  coord_sf(ylim = c(-58, 33), xlim = c(-138, -8)) +
  theme_void() + 
  theme(plot.background = element_rect(fill = "transparent", color = NA))

bio_NT_map

bio_AT_map <- ggplot() + 
  geom_sf(data = bio_AT, fill = "lightgrey", color = "lightgrey") +
  coord_sf(ylim = c(-40.5, 26.5), xlim = c(-15, 77.6)) +  
  theme_void() + 
  theme(plot.background = element_rect(fill = "transparent", color = NA))

bio_AT_map

bio_EP_map <- ggplot() + 
  geom_sf(data = bio_EP, fill = "lightgrey", color = "lightgrey") +
  coord_sf(ylim = c(24.5, 80), xlim = c(45, 174)) +  
  theme_void() + 
  theme(plot.background = element_rect(fill = "transparent", color = NA))

bio_EP_map

bio_WP_map <- ggplot() + 
  geom_sf(data = bio_WP, fill = "lightgrey", color = "lightgrey") +
  coord_sf(ylim = c(16, 80), xlim = c(-47, 90)) +  
  theme_void() + 
  theme(plot.background = element_rect(fill = "transparent", color = NA))

bio_WP_map

bio_AA_map <- ggplot() + 
  geom_sf(data = bio_AA, fill = "lightgrey", color = "lightgrey") +
  coord_sf(ylim = c(-54.75, 4.8), xlim = c(83, 180)) +  
  theme_classic() +
  theme_void() + 
  theme(plot.background = element_rect(fill = "transparent", color = NA))

bio_AA_map

bio_IM_map <- ggplot() + 
  geom_sf(data = bio_IM, fill = "lightgrey", color = "lightgrey") +
  coord_sf(ylim = c(-12.5, 36), xlim = c(65.85, 134)) +  
  theme_void() + 
  theme(plot.background = element_rect(fill = "transparent", color = NA))

bio_IM_map

bio_realms_map <- ggplot() + 
  geom_sf(data = realms, fill = "lightgrey", color = "lightgrey") +
  coord_sf(ylim = c(-56, 84), xlim = c(-180, 180)) +  
  theme_void() + 
  theme(plot.background = element_rect(fill = "transparent", color = NA))

bio_realms_map

make_hist_log <- function(data, fill_color = "#E31A1C", title = "", binwidth = 0.1) {
  data <- data[data$n > 0, ]
  data$log_n <- log(data$n)
  
  ggplot(data, aes(x = log_n)) +
    geom_histogram(fill = fill_color, color = fill_color, linewidth = 0.5, binwidth = binwidth) +
    coord_cartesian(xlim = c(0, 12)) +
    scale_x_continuous(
      breaks = seq(0, 12, by = 2),  # Increments of 2
      labels = scales::number_format(accuracy = 1)
    ) +
    scale_y_continuous(
      breaks = scales::pretty_breaks()) +
    labs(
      title = title,
      x = expression(log*"(Number of records)"),
      y = "Frequency"
    ) +
    theme_test(base_size = 7) +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 1),
      axis.text = element_text(size = 7),
      axis.title = element_text(size = 7)
    )
}

hist_NA <- make_hist_log(NA_counts, fill_color = "#A6CEE3") + ylab("") + xlab ("")
hist_NT <- make_hist_log(NT_counts, fill_color = "#FDBF6F") + ylab ("Frequency") + xlab ("Number of records (log)")
hist_WP <- make_hist_log(WP_counts, fill_color = "#FB9A99") + ylab ("") + xlab ("")
hist_EP <- make_hist_log(EP_counts, fill_color = "#B2DF8A") + ylab ("") + xlab ("Number of records (log)")
hist_IM <- make_hist_log(IM_counts, fill_color = "#FF7F00") + ylab ("") + xlab ("Number of records (log)")
hist_AT <- make_hist_log(AT_counts, fill_color = "#1F78B4") + ylab ("") + xlab ("Number of records (log)")
hist_AA <- make_hist_log(AA_counts, fill_color = "#E31A1C") + ylab ("") + xlab ("")
hist_GBIF <- make_hist_log(GBIF_counts, fill_color = "#999999") + ylab ("Frequency") + xlab ("")


# Plot histograms together with domain silhouettes

NA_plot <- ggdraw() +
  draw_plot(hist_NA) +
  draw_plot(bio_NA_map, x = 0.7, y = 0.6, width = 0.3, height = 0.3)
NA_plot

NT_plot <- ggdraw() +
  draw_plot(hist_NT) +
  draw_plot(bio_NT_map, x = 0.7, y = 0.6, width = 0.3, height = 0.3)
NT_plot

WP_plot <- ggdraw() +
  draw_plot(hist_WP) +
  draw_plot(bio_WP_map, x = 0.7, y = 0.6, width = 0.3, height = 0.3)
WP_plot

EP_plot <- ggdraw() +
  draw_plot(hist_EP) +
  draw_plot(bio_EP_map, x = 0.65, y = 0.6, width = 0.3, height = 0.3)
EP_plot

IM_plot <- ggdraw() +
  draw_plot(hist_IM) +
  draw_plot(bio_IM_map, x = 0.7, y = 0.6, width = 0.3, height = 0.3)
IM_plot

AT_plot <- ggdraw() +
  draw_plot(hist_AT) +
  draw_plot(bio_AT_map, x = 0.7, y = 0.6, width = 0.3, height = 0.3)
AT_plot

AA_plot <- ggdraw() +
  draw_plot(hist_AA) +
  draw_plot(bio_AA_map, x = 0.7, y = 0.6, width = 0.3, height = 0.3)
AA_plot

GBIF_plot <- ggdraw() +
  draw_plot(hist_GBIF) +
  draw_plot(bio_realms_map, x = 0.58, y = 0.57, width = 0.4, height = 0.4)
GBIF_plot

hists <- ggarrange(GBIF_plot, NA_plot, WP_plot, AA_plot, NT_plot, AT_plot, EP_plot, IM_plot, 
                     ncol = 4, nrow = 2, labels = c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)", "(g)", "(h)"),
                     font.label = list(size = 7, face = "plain", color ="black"))
hists 


# Cumulative number  Figure 1 --------------------------------------------------

# The code below is used to calculate the number of occurrence records 
# collected each year and the cumulative number of observed species

# Remove records lacking collection year in each biogeographic realm

Records_AA_year <- Records_AA %>%
  filter(!is.na(year))

Records_AT_year <- Records_AT %>%
  filter(!is.na(year))

Records_EP_year <- Records_EP %>%
  filter(!is.na(year))

Records_IM_year <- Records_IM %>%
  filter(!is.na(year))

Records_NA_year <- Records_NA %>%
  filter(!is.na(year))

Records_NT_year <- Records_NT %>%
  filter(!is.na(year))

Records_OC_year <- Records_OC %>%
  filter(!is.na(year))

Records_WP_year <- Records_WP %>%
  filter(!is.na(year))

# Nearctic ---------------------------------------------------------------------

Records_NA_year$duplicated <- duplicated(Records_NA_year$accepted_name) # assign a false to each new species 
Records_NA_year$specacum <- ifelse(Records_NA_year$duplicated == FALSE, 1, 0)# assign "1" to each new species and 0 to duplicated
Records_NA_year$recs <- as.integer(1) # assign 1 to each occurrence record
Records_NA_year_year <- aggregate(recs ~ year, Records_NA_year, sum) # Sum species by year

Acum_NA_year <- subset(Records_NA_year[, c('accepted_name', 'year', 'specacum')], specacum == 1)
Acum_NA_year$Counts <- as.integer(1) # assign 1 to each occurrence record

Acum_NA_year_2 <- aggregate(Counts ~ year, Acum_NA_year, sum) # Sum species by year
Acum_NA_year_2 <- aggregate(cbind(Counts, specacum)  ~ year, Acum_NA_year, sum) # groups by year
Acum_NA_year_2$acumulado <- cumsum(Acum_NA_year_2$specacum) # sums the number of species described each year

Records_NA_year_year <- merge(Records_NA_year_year, Acum_NA_year_2, by = 'year', all.x = TRUE)

# Afrotropical -------------------------------------------------------------------

Records_AT_year$duplicated <- duplicated(Records_AT_year$accepted_name) # assign a false to each new species 
Records_AT_year$specacum <- ifelse(Records_AT_year$duplicated == FALSE, 1, 0)# assign "1" to each new species and 0 to duplicated
Records_AT_year$recs <- as.integer(1) # assign 1 to each occurrence record
Records_AT_year_year <- aggregate(recs ~ year, Records_AT_year, sum) # Sum species by year

Acum_AT_year <- subset(Records_AT_year[, c('accepted_name', 'year', 'specacum')], specacum == 1)
Acum_AT_year$Counts <- as.integer(1) # assign 1 to each occurrence record

Acum_AT_year_2 <- aggregate(Counts ~ year, Acum_AT_year, sum) # Sum species by year
Acum_AT_year_2 <- aggregate(cbind(Counts, specacum)  ~ year, Acum_AT_year, sum) # groups by year
Acum_AT_year_2$acumulado <- cumsum(Acum_AT_year_2$specacum) # sums the number of species described each year

Records_AT_year_year <- merge(Records_AT_year_year, Acum_AT_year_2, by = 'year', all.x = TRUE)

# Eastern Palearctic -----------------------------------------------------------

Records_EP_year$duplicated <- duplicated(Records_EP_year$accepted_name) # assign a false to each new species 
Records_EP_year$specacum <- ifelse(Records_EP_year$duplicated == FALSE, 1, 0)# assign "1" to each new species and 0 to duplicated
Records_EP_year$recs <- as.integer(1) # assign 1 to each occurrence record
Records_EP_year_year <- aggregate(recs ~ year, Records_EP_year, sum) # Sum species by year

Acum_EP_year <- subset(Records_EP_year[, c('accepted_name', 'year', 'specacum')], specacum == 1)
Acum_EP_year$Counts <- as.integer(1) # assign 1 to each occurrence record

Acum_EP_year_2 <- aggregate(Counts ~ year, Acum_EP_year, sum) # Sum species by year
Acum_EP_year_2 <- aggregate(cbind(Counts, specacum)  ~ year, Acum_EP_year, sum) # groups by year
Acum_EP_year_2$acumulado <- cumsum(Acum_EP_year_2$specacum) # sums the number of species described each year

Records_EP_year_year <- merge(Records_EP_year_year, Acum_EP_year_2, by = 'year', all.x = TRUE)

# Indomalayan ------------------------------------------------------------------

Records_IM_year$duplicated <- duplicated(Records_IM_year$accepted_name) # assign a false to each new species 
Records_IM_year$specacum <- ifelse(Records_IM_year$duplicated == FALSE, 1, 0)# assign "1" to each new species and 0 to duplicated
Records_IM_year$recs <- as.integer(1) # assign 1 to each occurrence record
Records_IM_year_year <- aggregate(recs ~ year, Records_IM_year, sum) # Sum species by year

Acum_IM_year <- subset(Records_IM_year[, c('accepted_name', 'year', 'specacum')], specacum == 1)
Acum_IM_year$Counts <- as.integer(1) # assign 1 to each occurrence record

Acum_IM_year_2 <- aggregate(Counts ~ year, Acum_IM_year, sum) # Sum species by year
Acum_IM_year_2 <- aggregate(cbind(Counts, specacum)  ~ year, Acum_IM_year, sum) # groups by year
Acum_IM_year_2$acumulado <- cumsum(Acum_IM_year_2$specacum) # sums the number of species described each year

Records_IM_year_year <- merge(Records_IM_year_year, Acum_IM_year_2, by = 'year', all.x = TRUE)

# Australasia ------------------------------------------------------------------

Records_AA_year$duplicated <- duplicated(Records_AA_year$accepted_name) # assign a false to each new species 
Records_AA_year$specacum <- ifelse(Records_AA_year$duplicated == FALSE, 1, 0)# assign "1" to each new species and 0 to duplicated
Records_AA_year$recs <- as.integer(1) # assign 1 to each occurrence record
Records_AA_year_year <- aggregate(recs ~ year, Records_AA_year, sum) # Sum species by year

Acum_AA_year <- subset(Records_AA_year[, c('accepted_name', 'year', 'specacum')], specacum == 1)
Acum_AA_year$Counts <- as.integer(1) # assign 1 to each occurrence record

Acum_AA_year_2 <- aggregate(Counts ~ year, Acum_AA_year, sum) # Sum species by year
Acum_AA_year_2 <- aggregate(cbind(Counts, specacum)  ~ year, Acum_AA_year, sum) # groups by year
Acum_AA_year_2$acumulado <- cumsum(Acum_AA_year_2$specacum) # sums the number of species described each year

Records_AA_year_year <- merge(Records_AA_year_year, Acum_AA_year_2, by = 'year', all.x = TRUE)

# Neotropics --------------------------------------------------------------------

Records_NT_year$duplicated <- duplicated(Records_NT_year$accepted_name) # assign a false to each new species 
Records_NT_year$specacum <- ifelse(Records_NT_year$duplicated == FALSE, 1, 0)# assign "1" to each new species and 0 to duplicated
Records_NT_year$recs <- as.integer(1) # assign 1 to each occurrence record
Records_NT_year_year <- aggregate(recs ~ year, Records_NT_year, sum) # Sum species by year

Acum_NT_year <- subset(Records_NT_year[, c('accepted_name', 'year', 'specacum')], specacum == 1)
Acum_NT_year$Counts <- as.integer(1) # assign 1 to each occurrence record

Acum_NT_year_2 <- aggregate(Counts ~ year, Acum_NT_year, sum) # Sum species by year
Acum_NT_year_2 <- aggregate(cbind(Counts, specacum)  ~ year, Acum_NT_year, sum) # groups by year
Acum_NT_year_2$acumulado <- cumsum(Acum_NT_year_2$specacum) # sums the number of species described each year

Records_NT_year_year <- merge(Records_NT_year_year, Acum_NT_year_2, by = 'year', all.x = TRUE)

# Western Palearctic -----------------------------------------------------------

Records_WP_year$duplicated <- duplicated(Records_WP_year$accepted_name) # assign a false to each new species 
Records_WP_year$specacum <- ifelse(Records_WP_year$duplicated == FALSE, 1, 0)# assign "1" to each new species and 0 to duplicated
Records_WP_year$recs <- as.integer(1) # assign 1 to each occurrence record
Records_WP_year_year <- aggregate(recs ~ year, Records_WP_year, sum) # Sum species by year

Acum_WP_year <- subset(Records_WP_year[, c('accepted_name', 'year', 'specacum')], specacum == 1)
Acum_WP_year$Counts <- as.integer(1) # assign 1 to each occurrence record

Acum_WP_year_2 <- aggregate(Counts ~ year, Acum_WP_year, sum) # Sum species by year
Acum_WP_year_2 <- aggregate(cbind(Counts, specacum)  ~ year, Acum_WP_year, sum) # groups by year
Acum_WP_year_2$acumulado <- cumsum(Acum_WP_year_2$specacum) # sums the number of species described each year

Records_WP_year_year <- merge(Records_WP_year_year, Acum_WP_year_2, by = 'year', all.x = TRUE)

# Oceania ----------------------------------------------------------------------

Records_OC_year$duplicated <- duplicated(Records_OC_year$accepted_name) # assign a false to each new species 
Records_OC_year$specacum <- ifelse(Records_OC_year$duplicated == FALSE, 1, 0)# assign "1" to each new species and 0 to duplicated
Records_OC_year$recs <- as.integer(1) # assign 1 to each occurrence record
Records_OC_year_year <- aggregate(recs ~ year, Records_OC_year, sum) # Sum species by year

Acum_OC_year <- subset(Records_OC_year[, c('accepted_name', 'year', 'specacum')], specacum == 1)
Acum_OC_year$Counts <- as.integer(1) # assign 1 to each occurrence record

Acum_OC_year_2 <- aggregate(Counts ~ year, Acum_OC_year, sum) # Sum species by year
Acum_OC_year_2 <- aggregate(cbind(Counts, specacum)  ~ year, Acum_OC_year, sum) # groups by year
Acum_OC_year_2$acumulado <- cumsum(Acum_OC_year_2$specacum) # sums the number of species described each year

Records_OC_year_year <- merge(Records_OC_year_year, Acum_OC_year_2, by = 'year', all.x = TRUE)

# Figure ------------------------------------------------------------------

plotyDouble <- function(data, x, y, z, x_breaks = NULL, y_breaks_by = 20, 
                        x_lab = "Year", y_lab = "Number of records", y2_lab = "OBserved number of species") {
  ggplot() + 
    geom_bar(mapping = aes(x = x, y = y * 1), stat = "identity", fill = 'black') +
    geom_line(data, 
              mapping = aes(x = year, y = acumulado * prop), 
              linewidth = 0.5, color = 'darkgrey') +
    
    scale_x_continuous(name = x_lab, 
                       breaks = if (!is.null(x_breaks)) x_breaks else waiver()) +
    
    scale_y_continuous(name = y_lab,
                       sec.axis = sec_axis(~ ./prop, name = y2_lab, 
                                           breaks = seq(0, max(y), by = y_breaks_by))) +
    
    theme_test(base_size = 6) +
    theme(axis.line = element_line(colour = "black", linewidth = 0.1, linetype = "solid"), 
          axis.title.x = element_text(vjust = -0.5, size = 6, color = "black"),
          axis.title.y = element_text(vjust = 2, size = 6),
          axis.title = element_text(size = 6),
          axis.text = element_text(size = 6, color = "black"), 
          axis.ticks.length = unit(0.10, "cm"))
}

sacs_AA <- Records_AA_year_year[!is.na(Records_AA_year_year$acumulado),]
prop <- round(max(Records_AA_year_year$recs)/max(sacs_AA$acumulado), -1)
sacRecs_AA <- plotyDouble(sacs_AA, 
                          x = Records_AA_year_year$year, 
                          y = Records_AA_year_year$recs, 
                          z = Records_AA_year_year$acumulado,
                          x_breaks = seq(1770, 2025, by = 60),
                          y_breaks_by = 60,
                          y_lab = "",
                          x_lab = "",
                          y2_lab = "")

AA_sac <- ggdraw() +
  draw_plot(sacRecs_AA) +
  draw_plot(bio_AA_map, x = 0.2, y = 0.65, width = 0.3, height = 0.3)
AA_sac

sacs_AT <- Records_AT_year_year[!is.na(Records_AT_year_year$acumulado),]
prop <- round(max(Records_AT_year_year$recs)/max(sacs_AT$acumulado), -1)
sacRecs_AT <- plotyDouble(sacs_AT, 
                          x = Records_AT_year_year$year, 
                          y = Records_AT_year_year$recs, 
                          z = Records_AT_year_year$acumulado,
                          x_breaks = seq(1850, 2020, by = 40),
                          y_breaks_by = 90,
                          y_lab = "Number of records",
                          x_lab = "Year",
                          y2_lab = "")
AT_sac <- ggdraw() +
  draw_plot(sacRecs_AT) +
  draw_plot(bio_AT_map, x = 0.25, y = 0.65, width = 0.3, height = 0.3)
AT_sac

sacs_NT <- Records_NT_year_year[!is.na(Records_NT_year_year$acumulado),]
prop <- round(max(Records_NT_year_year$recs)/max(sacs_NT$acumulado), -1)
sacRecs_NT <- plotyDouble(sacs_NT, 
                          x = Records_NT_year_year$year, 
                          y = Records_NT_year_year$recs, 
                          z = Records_NT_year_year$acumulado,
                          x_breaks = seq(1800, 2025, by = 75),
                          y_breaks_by = 230,
                          y_lab = "",
                          x_lab = "",
                          y2_lab = "Observed number of species")
NT_sac <- ggdraw() +
  draw_plot(sacRecs_NT) +
  draw_plot(bio_NT_map, x = 0.2, y = 0.65, width = 0.3, height = 0.3)
NT_sac


sacs_IM <- Records_IM_year_year[!is.na(Records_IM_year_year$acumulado),]
prop <- round(max(Records_IM_year_year$recs)/max(sacs_IM$acumulado), -1)
sacRecs_IM <- plotyDouble(sacs_IM, 
                          x = Records_IM_year_year$year, 
                          y = Records_IM_year_year$recs, 
                          z = Records_IM_year_year$acumulado,
                          x_breaks = seq(1845, 2025, by = 45),
                          y_breaks_by = 130,
                          y_lab = "",
                          x_lab = "Year",
                          y2_lab = "Observed number of species")
IM_sac <- ggdraw() +
  draw_plot(sacRecs_IM) +
  draw_plot(bio_IM_map, x = 0.25, y = 0.65, width = 0.3, height = 0.3)
IM_sac

sacs_EP <- Records_EP_year_year[!is.na(Records_EP_year_year$acumulado),]
prop <- round(max(Records_EP_year_year$recs)/max(sacs_EP$acumulado), -1)
sacRecs_EP <- plotyDouble(sacs_EP, 
                          x = Records_EP_year_year$year, 
                          y = Records_EP_year_year$recs, 
                          z = Records_EP_year_year$acumulado,
                          x_breaks = seq(1882, 2022, by = 35),
                          y_breaks_by = 35,
                          y_lab = "",
                          x_lab = "Year",
                          y2_lab = "")
EP_sac <- ggdraw() +
  draw_plot(sacRecs_EP) +
  draw_plot(bio_EP_map, x = 0.25, y = 0.65, width = 0.3, height = 0.3)
EP_sac

sacs_NA <- Records_NA_year_year[!is.na(Records_NA_year_year$acumulado),]
prop <- round(max(Records_NA_year_year$recs)/max(sacs_NA$acumulado), -1)
sacRecs_NA <- plotyDouble(sacs_NA, 
                          x = Records_NA_year_year$year, 
                          y = Records_NA_year_year$recs, 
                          z = Records_NA_year_year$acumulado,
                          x_breaks = seq(1700, 2020, by = 80),
                          y_breaks_by = 50,
                          y_lab = "Number of records",
                          x_lab = "",
                          y2_lab = "")
NA_sac <- ggdraw() +
  draw_plot(sacRecs_NA) +
  draw_plot(bio_NA_map, x = 0.3, y = 0.65, width = 0.3, height = 0.3)
NA_sac


sacs_WP <- Records_WP_year_year[!is.na(Records_WP_year_year$acumulado),]
prop <- round(max(Records_WP_year_year$recs)/max(sacs_WP$acumulado), -1)
sacRecs_WP <- plotyDouble(sacs_WP, 
                          x = Records_WP_year_year$year, 
                          y = Records_WP_year_year$recs, 
                          z = Records_WP_year_year$acumulado,
                          x_breaks = seq(1575, 2025, by = 110),
                          y_breaks_by = 22,
                          y_lab = "",
                          x_lab = "",
                          y2_lab = "")
WP_sac <- ggdraw() +
  draw_plot(sacRecs_WP) +
  draw_plot(bio_WP_map, x = 0.2, y = 0.65, width = 0.3, height = 0.3)
WP_sac


acc <- ggarrange(NA_sac, WP_sac, AA_sac, NT_sac, AT_sac, EP_sac, IM_sac,
                     ncol = 4, nrow = 2, labels = c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)", "(g)"),
                     font.label = list(size = 6.5, face = "plain", color ="black"), 
                     widths = rep(1,4), heights = rep(1, 2),  hjust = -0.1,  
                     vjust = 1.5, align = "hv")
acc 






