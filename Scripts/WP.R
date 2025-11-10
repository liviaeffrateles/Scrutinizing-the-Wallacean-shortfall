library(rgdal)  
library(sp)   
library(raster)   
library(maps)
library(sf)
library(nFactors)
library(psych)
library(colorRamps)  
library(ggplot2)
library(mapview)
library(terra)
library(rnaturalearth)
library(dplyr)
library(openxlsx)
library(stringr)

# Import the 19 CHELSA bioclimatic variables (1° resolution) cropped to the Afrotropical realm

#Raster
path_env <- "D:/Scrutinizing-the-Wallacean-shortfall/Chelsa"
files_env <- list.files(path_env, pattern = "^climate_W_Palearctic.*\\.tif$", full.names = TRUE)
allrasters <- lapply(files_env, raster)
clim_env <- raster::stack(files_env)

# Check the resulting aggregated raster to confirm the new resolution
res(clim_env)
r_WPA <- raster(clim_env, 1) #extract 1 raster to check NEW resolution of cells
res_WPA <- res(r_WPA)  # Save value for later analysis
climate_W_Palearctic <- clim_env # create climate stack of our new rasters
plot(climate_W_Palearctic[[1]])

# PCA ---------------------------------------------------------------------

library(psych)
library(maps)
library(colorRamps)

# We reduce all the variables to fewer variables using a PCA.
# 1) Here we standardize and prepare the data.
valuesW_Palearctic <- as.data.frame(values(climate_W_Palearctic)) # get env values from rasters

# 2) Standardization of the variables (Normalization - Mean = 0 and Std = 1)
std <- function(x){(x - mean(x, na.rm = T)) / sd(x, na.rm = T)} 
valuesWPalearcticStd <- apply(valuesW_Palearctic, 2, std) # apply to dataframe of climate values

rem <- apply(is.na(valuesWPalearcticStd), 1, any)
PCAWPalearctic <- as.data.frame(valuesWPalearcticStd[!rem, ])

# 3) Run the PCA (2 axis)
matWPalearctic <- matrix(runif(nrow(PCAWPalearctic) * ncol(PCAWPalearctic), 0.00001, 0.00009), 
                         ncol = ncol(PCAWPalearctic)) # add a very small randomness to avoid singularity
PCAWPalearctic2 <- PCAWPalearctic + matWPalearctic
myPCAWPalearctic <- principal(PCAWPalearctic2,
                              nfactors = 2,
                              rotate = "varimax",
                              scores = T)

# 4)
loadings_WP <- as.data.frame(unclass(myPCAWPalearctic$loadings))

loadings_PC1_WP <- loadings_WP %>%
  mutate(abs_RC1 = abs(RC1)) %>%
  arrange(desc(abs_RC1))

loadings_PC2_WP <- loadings_WP %>%
  mutate(abs_RC2 = abs(RC2)) %>%
  arrange(desc(abs_RC2))

# 4) Plot PCA
biplot.psych(myPCAWPalearctic, xlim.s = c(-6, 2), ylim.s = c(-3, 3)) 
rownames(loadings_WP) <- str_extract(rownames(loadings_WP), "bio\\d+")
loadings_WP <- tibble::rownames_to_column(loadings_WP, var = "variables")
scores_WP <- myPCAWPalearctic$scores
scale_factor <- 2.5

biplot_WP <- ggplot() +
  geom_point(data = as.data.frame(scores_WP), aes(x = RC1, y = RC2), color = "black", size = 1) +
  geom_text(data = loadings_WP, aes(x = RC1  * scale_factor * 1.05, y = RC2  * scale_factor * 1.05, label = variables), 
            vjust = -0.5, hjust = 0.5, size = 2, color = "#F0A202", fontface = "bold") +
  geom_segment(data = loadings_WP, aes(x = 0, y = 0, xend = RC1 * scale_factor, yend = RC2 * scale_factor), 
               arrow = arrow(type = "open", length = unit(0.1, "inches")), color = "#0072B2") +
  scale_x_continuous(breaks = seq(-6, 2, by = 1)) + 
  scale_y_continuous(breaks = seq(-3, 3, by = 1)) +
  xlab(paste("PC1 (55%)")) +
  ylab(paste("PC2 (28%)")) +
  theme_test(base_size = 12) +
  theme(aspect.ratio = 1)

biplot_WP

# PCA of our study area ---------------------------------------------------

# Instead of assuming the absolute values of all Chelsa variables, 
# the next map assumes values of the linear combinations between these variables (PCA scores).
valuesWPalearctic_2 <- valuesWPalearcticStd[, 1:2] # create a vector with the same length as v2 but 2 columns
valuesWPalearctic_2[!rem, ] <- myPCAWPalearctic$scores 
climate_PCA_WPalearctic <- subset(climate_W_Palearctic, 1:2)
values(climate_PCA_WPalearctic) <- valuesWPalearctic_2 # insert pca values of v3 into raster
names(climate_PCA_WPalearctic) <- c("PC1", "PC2") # rename
# coords_PCA_WP <- as.data.frame(climate_PCA_WPalearctic, xy = T)
# write.table(coords_PCA_WP, "Processed_Data/coords_PCA_WP.txt") # Save for environmental distance calculation

# Plot study area by each axis from PCA values
world_sf <- ne_countries(scale = "medium", returnclass = "sf") %>%
  dplyr::filter(name != "Antarctica")

par(mfrow = c(1, 2))
plot(climate_PCA_WPalearctic, 1, font = 1, font.lab = 1, add = T)  # PCA axis 1
plot(world_sf, add = TRUE, border = "gray30", col = NA, lwd = 0.8, xlim = c(-47, 90), ylim = c(16, 80))
plot(climate_PCA_WPalearctic, 2, font = 1, font.lab = 1)  # PCA axis 2
plot(world_sf, add = TRUE, border = "gray30", col = NA, lwd = 0.8, xlim = c(-47, 90), ylim = c(16, 80))

# Environmental space -----------------------------------------------------

# The next step is to create the environmental space using the two PCA scores.
# Create a Cartesian plan with PC1 and PC2 scores

climate_PCA_values_WPalearctic <- values(climate_PCA_WPalearctic) # get env values from PCA

# Transform the two vars we want into the env space, 
# by taking the min and max scores (PCA score values) for each PCA axis
# and create a raster object. 
xminWPA <- min(climate_PCA_values_WPalearctic[, 1], na.rm = TRUE) 
xmaxWPA <- max(climate_PCA_values_WPalearctic[, 1], na.rm = TRUE)
yminWPA <- min(climate_PCA_values_WPalearctic[, 2], na.rm = TRUE)
ymaxWPA <- max(climate_PCA_values_WPalearctic[, 2], na.rm = TRUE)

# This function creates the cartesian plan comprising the min and max PCA scores

env_space_WPA <- raster(xmn = xminWPA, xmx = xmaxWPA, 
                        ymn = yminWPA, ymx = ymaxWPA, crs = 4326,
                        res = 0.5) # Resolution of the env cell can change
ncell(env_space_WPA)
values(env_space_WPA) <- 0
env_space_area_WPA <- env_space_WPA # duplicate this object for the next step

### Insert our PCA values into this env_space
# Extract PC1 and PC2, convert it in "classes of values" to plot in the map
env_space_WPA_v <- raster::extract(env_space_WPA, 
                                   y = na.omit(climate_PCA_values_WPalearctic), 
                                   cellnumbers = TRUE)[, 1]
table(env_space_WPA_v)
n_env_space_WPA_v <- table(env_space_WPA_v) # this function counts the frequency of "climates" (PC scores)

values(env_space_area_WPA)[as.numeric(names(n_env_space_WPA_v))] <- n_env_space_WPA_v
area_values_WPA <- values(env_space_area_WPA)
area_values_WPA[area_values_WPA == 0] <- NA
values(env_space_area_WPA) <- area_values_WPA
pol <- rasterToPolygons(env_space_area_WPA, dissolve = FALSE)

# The following map shows the frequency of climates in the study area.
# The axes shows all the PC1 and PC2 scores. 
# Each cell "roughly" represents a climate type.
# The colors indicate the frequency of climate types.

cols <- colorRampPalette(c("#56B4E9", "#33A02C", "#FDE725FF", "#CC79A7", "#D55E00"))

plot(env_space_area_WPA, 
     xlab = "PC1", 
     ylab = "PC2",
     col = cols(1000),
     xlim = c(-5.5, 2.5),
     ylim = c(-3.5, 3.5),
     font = 2, font.lab = 2,
     cex.axis = 0.8,
     cex.lab = 0.8,
     axis.args = list(cex.axis = 0.8)
)
plot(pol, add = TRUE)

# writeRaster(env_space_area_WPA, "Processed_Data/env_space_WP.tif", format = "GTiff", overwrite = TRUE)

# Well-sampled cells ------------------------------------------------------

# Import well-sampled cells with a slope threshold of 0.05 from the inventory completeness analysis
# well_sampled_1 <- read.csv("Data/well_sampled_1.csv")

# Select only well-sampled cells within the boundaries of Western Palearctic realm
well_sampled_WPA_1 <- well_sampled_1 %>% 
  dplyr::filter(Realm == "WP")

#_______________________________________________________________________________
# Now we calculate the environmental space for our Well-sampled cells
# Using estimators output from 'InventoryCompleteness' script:
well_sampled_WPA_1 <- well_sampled_WPA_1[, c('Longitude','Latitude')]
# Plot the well-sampled occurrences on the map
raster::plot(r_WPA, main = "", col = "gray", legend = FALSE)
points(well_sampled_WPA_1, col = rgb(1, 0, 0, .5), pch = 20, cex = 2)

# Extract climate values of Well-sampled cells
values_WS_WPA <- raster::extract(climate_W_Palearctic, 
                                 well_sampled_WPA_1,
                                 cellnumbers = TRUE)[, 1]

# write.table(values_WS_WPA, "Processed_Data/values_WS_WP.txt") # Save for environmental distance calculation

coords_WS_WPA <- climate_PCA_values_WPalearctic[values_WS_WPA, ]

# write.table(coords_WS_WPA, "Processed_Data/coords_WS_WP.txt") # Save for environmental distance calculation
cell_WS_WPA <- raster::extract(env_space_WPA, 
                               coords_WS_WPA,
                               cellnumbers = TRUE)[, 1]

n_WS_WPA <- table(cell_WS_WPA)
env_space_WS_WPA <- env_space_WPA
values(env_space_WS_WPA)[as.numeric(names(n_WS_WPA))] <- n_WS_WPA
ncell(cell_WS_WPA)

# The following map shows the frequency of climates in the well-sampled cells 
# in the study area (legend of this figure is continuous BUT values are integer numbers, 
# so manually adapt the legend to each case)

{tiff("Figures/Climate_Types_WS.tif", width = 10, height = 10, res = 600, units = "cm")
plot(env_space_WS_WPA, 
     xlab = "PC1", 
     ylab = "PC2",
     cex.lab = 1,
     col = c("transparent", cols(1000)),
     xlim = c(-5.5, 1.7),
     ylim = c(-2.9, 2.7),
     font = 2, font.lab = 2,
     cex.axis = 0.8,
     cex.lab = 0.8,
     axis.args = list(cex.axis = 0.8))
plot(pol, add = TRUE)
}
dev.off()

# writeRaster(env_space_WS_WPA, "Processed_Data/env_space_WS_WP.tif", format = "GTiff", overwrite = TRUE)

# Schoener's D ------------------------------------------------------------

# Schoener's D: quantifies the overlap between the location of well-sampled 
# sites and cells with most frequent conditions 
# The Schoener's D index varies from zero (total lack of congruence) to one 
# (total congruence) 

# Transform the abundance of each cell into probabilities.
area_values_WPA <- area_values_WPA/sum(area_values_WPA, na.rm = TRUE) # Relative frequency of climate type for all the study area
WS_values_WPA <- values(env_space_WS_WPA)
WS_values_WPA <- WS_values_WPA/sum(WS_values_WPA, na.rm = TRUE) # Relative frequency of climate type for well-sampled cells

# Calculate the climate overlap using Schoener's D.
# D values close to 1 indicate that the location of well-sampled 
# sites coincide with climate conditions 
# that are frequently found in the study area
SchoenersD <- function(x, y) {
  sub_values <- abs(x - y)
  D <- 1 - (sum(sub_values, na.rm = TRUE) / 2)
  return(D)
}

D_WPA <- SchoenersD(area_values_WPA, WS_values_WPA)
print(paste("Climate overlap between well-sampled cells and the study area, given by the observed Schoener's D equals = ", 
            round(D_WPA, 3), "%"))

# We create a null model to test if D values is different from a 
# random distribution of D values calculated from randomly sampling occurrence 
# records.
set.seed(0)
replications <- 1000 # Choose the number of replications
D_rnd <- numeric(replications)
for (i in 1:replications) {
  rnd <- sample(env_space_WPA_v, length(well_sampled_WPA_1), replace = TRUE)
  n_rnd <- table(rnd)
  env_space_rnd <- env_space_WPA
  values(env_space_rnd)[as.numeric(names(n_rnd))] <- n_rnd
  rnd_values <- values(env_space_rnd)
  rnd_values <- rnd_values/sum(rnd_values, na.rm = TRUE)
  D_rnd[i] <- SchoenersD(area_values_WPA, rnd_values)
}

# p-value from previous analysis
# If p < 0.05, it means that the location of well-sampled sites does not 
# coincide with areas with climate conditions frequently found in your study area
p_WPA <- (sum(D_WPA > D_rnd) + 1) / (length(D_rnd) + 1) # Unicaudal test
print(paste("p value equals = ", round(p_WPA, 3)))

# Plot distribution of model values
hist(D_rnd, 10, # Write number of bins
     xlim = c(0, 1), 
     main = (""),
     xlab = "D", 
     col = rgb(.5, .5, .5), 
     border = FALSE,
     font = 1, font.lab = 1)
abline(v = D_WPA, col = "red", lty = 2) # Red line show the observed Schoeners' D value

############ Kruskal-Wallis test ##################################
# The following map shows the distribution of well-sampled cells (red bars) 
# vs the study area cells (gray bars)
### For a) PC1
  hist(myPCAWPalearctic$scores[, 1], 
       breaks = ncol(env_space_WPA),
       freq = F,
       col = "grey", 
       border = FALSE, 
       xlim= c(-6, 2.03), 
       ylim = c(0, 0.8), 
       main = "",
       xlab = "PC1", 
       font = 2, font.lab = 2, 
       cex.lab = 0.8, cex.axis = 0.8)
  lines(density(na.omit(myPCAWPalearctic$scores[, 1])), col = "black", lwd = 2)
  
  hist(coords_WS_WPA[, 1],
       add = TRUE,
       col = adjustcolor("#0072B2", alpha.f = 0.5),
       freq = F,
       border = FALSE)
  lines(density(na.omit(coords_WS_WPA[, 1])), col = "#0072B2", lwd = 2)

# X axis is a probability density
# Kruskal-Wallis verifies whether 1) the distribution of well-sampled sites 
# is an unbiased subset of the entire climate conditions of the Atlantic forest. 
# If this is so, p > 0.05 
x_1_WPA <- c(myPCAWPalearctic$scores[, 1], coords_WS_WPA[, 1])
g_1_WPA <- as.factor(c(rep("area", length(myPCAWPalearctic$scores[, 1])),
                       rep("WS", length(coords_WS_WPA[, 1]))))
kruskal.test(x_1_WPA ~ g_1_WPA)
# Kolmogorov smirnov test###
ks.test(myPCAWPalearctic$scores[, 1], coords_WS_WPA[, 1])

### For b) PC2
  hist(myPCAWPalearctic$scores[, 2], 
       breaks = ncol(env_space_WPA),
       freq = F,
       col = "grey", 
       border = FALSE, 
       xlim = c(-2.8, 3.2), 
       ylim = c(0, 0.5), 
       main = "",
       xlab = "PC2", 
       font = 2, font.lab = 2, 
       cex.lab = 0.8, cex.axis = 0.8)
  lines(density(na.omit(myPCAWPalearctic$scores[, 2])), col = "black", lwd = 2)
  
  hist(coords_WS_WPA[, 2],
       breaks = ncol(env_space_WPA), 
       add = TRUE,
       col = adjustcolor("#0072B2", alpha.f = 0.5),
       freq = F,
       border = FALSE)
  lines(density(na.omit(coords_WS_WPA[, 2])), col = "#0072B2", lwd = 2)

x_2_WPA <- c(myPCAWPalearctic$scores[, 2], coords_WS_WPA[, 2])
g_2_WPA <- as.factor(c(rep("area", length(myPCAWPalearctic$scores[, 2])),
                       rep("WS", length(coords_WS_WPA[, 2]))))

kruskal.test(x_2_WPA ~ g_2_WPA)

# Kolmogorov smirnov test###
ks.test(myPCAWPalearctic$scores[, 2], coords_WS_WPA[, 2])

# Rarity ------------------------------------------------------------------

# We can check how many environmental space has been sampled 
# and how does these cells look like.
surface_WPA <- WS_values_WPA
surface_WPA[is.na(area_values_WPA) | area_values_WPA == 0] <- NA
sampled_WPA <- sum(surface_WPA > 0, na.rm = TRUE)/ sum(surface_WPA >= 0, na.rm = TRUE) * 100
percen_WPA <- round(sampled_WPA, 2)
prin_WPA <- paste0(percen_WPA, "% of our study area climate types covered by well-sampled cells")
print(prin_WPA)

# Is this env. space sampled corresponding to rare climates?
# First, we make values vary from 0 to 1 according to their rarity:
# Values close to 0, are very common, values close to 1 very rare
mini_WPA <- min(area_values_WPA, na.rm = TRUE) # less frequent value
# rarity index, also called Min-Max scalling
area_values01WPA <- abs(1 - (area_values_WPA - mini_WPA) / (max(area_values_WPA, na.rm = TRUE) - mini_WPA)) 
# see http://rasbt.github.io/mlxtend/user_guide/preprocessing/minmax_scaling/
hist(area_values01WPA, 
     breaks = ncol(env_space_WPA), 
     freq = FALSE, col = "white",
     main = "", 
     xlab = "Climate rarity", border = FALSE, 
     ylim = c(0, 10),
     font = 2, font.lab = 2,
     cex.lab = 0.8, cex.axis = 0.8)
lines(density(na.omit(area_values01WPA)), col = "black", lwd = 2)
lines(density(na.omit(area_values01WPA[surface_WPA > 0])), col = "#0072B2", lwd = 2)

# Kruskal-Wallis test to see if the distribution of rarities for the sampled 
# occurrences differs from the distribution observed in the entire area
x_WPA <- c(area_values01WPA, area_values01WPA[surface_WPA > 0])
g_WPA <- as.factor(c(rep("area", length(area_values01WPA)),
                     rep("ws", length(area_values01WPA[surface_WPA > 0]))))
kruskal.test(x_WPA ~ g_WPA)
# Kolmogorov smirnov test
area_values03_WPA <- (na.omit(area_values01WPA))
area_values04_WPA <- (na.omit(area_values01WPA[surface_WPA > 0]))
ks.test(area_values03_WPA, area_values04_WPA)

# This result indicates that the under-sampled area is composed 
# mainly by rare climates. You can inspect the environmental space figures 
# to check which areas were not sampled.

# Finally MAP the climatic rarity
rarity_env_WPA <- env_space_WPA
values(rarity_env_WPA) <- area_values01WPA
rarity_Percell_WPA <- raster::extract(rarity_env_WPA, climate_PCA_values_WPalearctic)
rarity_map_WPA <- r_WPA
values(rarity_map_WPA) <- rarity_Percell_WPA

cols <- colorRampPalette(c('#e0f3db','#a8ddb5','#43a2ca'))
x11()
plot(rarity_map_WPA, 
     col = cols(10), 
     font = 2, font.lab = 2, 
     ylim = c(14.91196, 81.90185), 
     xlim = c(-31.29, 69.07533)
)
points(well_sampled_WPA_1, col = "red", pch = 19, cex = 0.7)

# Convert raster to data frame
rarity_map_df_WPA <- as.data.frame(as(rarity_map_WPA, "SpatialPixelsDataFrame"))
colnames(rarity_map_df_WPA) <- c("value", "x", "y")

# write.csv(rarity_map_df_WPA, "Processed_Data/Rarity_WP.csv") # Save to create rarity map

