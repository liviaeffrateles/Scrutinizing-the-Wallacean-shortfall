library(colorRamps)
library(ggplot2)

# This script contains the code used to create the rarity figures and to depict well-sampled cells

# Import rarity data for each biogeographic realm
rarity_files <- list.files(path = "Processed_Data", pattern = "^Rarity_.*\\.csv$", full.names = TRUE)
rarity_list <- lapply(rarity_files, read.csv)
names(rarity_list) <- tools::file_path_sans_ext(basename(rarity_files))

# Import the well-sampled cells for NA, AA and WP
# WS_1 <- read.csv("Processed_Data/well_sampled_1.csv")

WS_1 <- WS_1 %>%
  mutate(Realm = ifelse(is.na(Realm), "NA", Realm))

WS_NA <- WS_1 %>%
  filter (Realm == "NA")

WS_AA <- WS_1 %>%
  filter (Realm == "AA")

WS_WP <- WS_1 %>%
  filter (Realm == "WP")

# Import the well-sampled cells for NT, AT, IM and EP
# WS_2 <- read.csv("Processed_Data/well_sampled_2.csv")

WS_IM <- WS_2 %>%
  filter (Realm == "IM")

WS_AT <- WS_2 %>%
  filter (Realm == "AT")

WS_NT <- WS_2 %>%
  filter (Realm == "NT")

WS_EP <- WS_2 %>%
  filter (Realm == "EP")

cols <- colorRampPalette(c("#9467BD", "#1F9E89", "#D8E219FF", "#FDE725"))

world_sf <- ne_countries(scale = "medium", returnclass = "sf") %>%
  dplyr::filter(name != "Antarctica")

# Rarity map for each biogeographic realm

create_rarity_map <- function(rarity_data, ws_data, xlim, ylim, legend_pos, legend_dir, guide_width, guide_height, legend_hjust = 0.5) {
  ggplot() +
    geom_tile(data = rarity_data, aes(x = x, y = y, fill = value)) +
    scale_fill_gradientn(colors = cols(10)) +
    geom_sf(data = world_sf, fill = NA, color = "black") +
    coord_sf(xlim = xlim, ylim = ylim) +
    geom_tile(data = ws_data, aes(x = Longitude, y = Latitude),
              width = 1, height = 1, linewidth = 0.2, color = "black", fill = NA) +
    theme_test(base_size = 10) +
    theme(
      panel.background = element_rect(fill = "#F7FBFF"),
      legend.background = element_rect(fill = NA, color = NA), 
      legend.key = element_rect(fill = NA, color = NA),
      plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "mm"),
      legend.text = element_text(size = 8), 
      legend.title = element_text(size = 8),
      axis.text = element_text(size = 10),
      axis.title.x = element_text(size = 10),
      axis.title.y = element_text(size = 10),
      axis.ticks.y = element_line(),
      legend.key.height = unit(0.4, 'cm'),
      legend.key.width = unit(0.4, 'cm'),
      legend.direction = legend_dir,
      legend.position = legend_pos
    ) +
    guides(fill = guide_colorbar(title.position = "top", title.hjust = legend_hjust,
                                 barwidth = guide_width, barheight = guide_height)) +
    labs(fill = "Rarity", x = "Longitude", y = "Latitude")
}


Rarity_IM_map <- create_rarity_map(rarity_list$Rarity_IM, WS_IM, c(64, 135), c(-12.5, 36),
                                   legend_pos = c(0.3, .13), legend_dir = "horizontal",
                                   guide_width = 5, guide_height = 0.5, legend_hjust = 0.5)

Rarity_AT_map <- create_rarity_map(rarity_list$Rarity_AT, WS_AT, c(-18, 78), c(-40.5, 25.9),
                                   legend_pos = c(0.15, 0.3), legend_dir = "vertical",
                                   guide_width = 0.5, guide_height = 4, legend_hjust = 0)

Rarity_NT_map <- create_rarity_map(rarity_list$Rarity_NT, WS_NT, c(-140, -7), c(-58.5, 31.9),
                                   legend_pos = c(0.1, .25), legend_dir = "vertical",
                                   guide_width = 0.5, guide_height = 4, legend_hjust = 0)

Rarity_EP_map <- create_rarity_map(rarity_list$Rarity_EP, WS_EP, c(45, 174), c(22, 79),
                                   legend_pos = c(0.8, .13), legend_dir = "horizontal",
                                   guide_width = 5, guide_height = 0.5, legend_hjust = 0.5)

Rarity_NA_map <- create_rarity_map(rarity_list$Rarity_NA, WS_NA, c(-167.2, -18), c(20, 83.6),
                                   legend_pos = c(0.1, .25), legend_dir = "vertical",
                                   guide_width = 0.5, guide_height = 4, legend_hjust = 0)

Rarity_AA_map <- create_rarity_map(rarity_list$Rarity_AA, WS_AA, c(75.7, 180), c(-58, 7.5),
                                   legend_pos = c(0.1, .25), legend_dir = "vertical",
                                   guide_width = 0.5, guide_height = 4, legend_hjust = 0)

Rarity_WP_map <- create_rarity_map(rarity_list$Rarity_WP, WS_WP, c(-47, 90), c(16, 80),
                                   legend_pos = c(0.1, .25), legend_dir = "vertical",
                                   guide_width = 0.5, guide_height = 4, legend_hjust = 0)

# Global Rarity map -------------------------------------------------------

bioRealms <- st_read("Grids/bioRealms.gpkg")

Rarity_map <- ggplot() +
    geom_sf(data = world_sf, fill = NA, color = "black") +
    geom_tile(data = rarity_list[["Rarity_NA"]], aes(x = x, y = y, fill = value)) + 
    geom_tile(data = rarity_list[["Rarity_IM"]], aes(x = x, y = y, fill = value)) +
    geom_tile(data = rarity_list[["Rarity_WP"]], aes(x = x, y = y, fill = value)) +
    geom_tile(data = rarity_list[["Rarity_EP"]], aes(x = x, y = y, fill = value)) +
    geom_tile(data = rarity_list[["Rarity_NT"]], aes(x = x, y = y, fill = value)) +
    geom_tile(data = rarity_list[["Rarity_AT"]], aes(x = x, y = y, fill = value)) +
    geom_tile(data = rarity_list[["Rarity_AA"]], aes(x = x, y = y, fill = value)) +
    scale_fill_gradientn(colors = cols(10)) +
    geom_sf(data = bioRealms, fill = NA, color = "#4D4D4D", size = 2) +
    geom_tile(data = WS_1, aes(x = Longitude, y = Latitude), width = 1, height = 1, linewidth = 0.12, color = "black", fill = NA) +
    geom_tile(data = WS_2, aes(x = Longitude, y = Latitude), width = 1, height = 1, linewidth = 0.12, color = "black", fill = NA) +
    coord_sf(xlim = c(-163.5, 163.5), ylim = c(-58, 75), expand = FALSE) +
    theme_test(base_size = 10) +
    theme(panel.background = element_rect(fill = "#F7FBFF"),
          legend.background = element_rect(fill = NA, color = NA), 
          legend.key = element_rect(fill = NA, color = NA))  +
    ggplot2::theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "mm"),
                   legend.text = element_text(size = 6.5),
                   legend.title = element_text(size = 6.5),
                   axis.text = element_text(size = 10),
                   axis.title.x = element_text(size = 10),
                   axis.title.y = element_text(size = 10),
                   axis.ticks.y = element_line(),
                   legend.key.height = unit(0.4, 'cm'),
                   legend.key.width = unit(0.4, 'cm'),
                   legend.position = c(0.1, .25),
                   legend.direction = "vertical") +
    guides(fill = guide_colorbar(title.position = "top", title.hjust = 0, 
                                 barwidth = 0.5, barheight = 4)) +
    labs(fill = "Rarity", x = "Longitude", y = "Latitude")



# Well-sampled cells  -----------------------------------------------------

  ws_cells <- ggplot() +
    geom_tile(data = WS_1, aes (x = Longitude, y = Latitude), width = 1, height = 1, linewidth = 0.2, color = "red", fill = NA) + 
    geom_tile(data = WS_2, aes (x = Longitude, y = Latitude), width = 1, height = 1, linewidth = 0.2, color = "red", fill = NA) +
    geom_sf(data = bioRealms, fill = NA, color = "#4D4D4D", size = 2) +
    coord_sf(xlim = c(-178.5, 178.5), ylim = c(-61.5, 81.7), expand = FALSE) +
    theme_test(base_size = 7) +
    theme(panel.background = element_rect(fill = "#F7FBFF")) +
    ggplot2::theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "mm"),
                   axis.text = element_text(size = 10),
                   axis.title.x = element_text(size = 10),
                   axis.title.y = element_text(size = 10),
                   axis.ticks.y = element_line()) +
    labs(x = "Longitude", y = "Latitude")

ws_cells


