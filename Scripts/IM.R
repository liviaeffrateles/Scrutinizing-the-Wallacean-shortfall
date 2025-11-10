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

# Import the 19 CHELSA bioclimatic variables (1° resolution) cropped to the Indomalayan realm

#Raster

path_env <- "D:/Scrutinizing-the-Wallacean-shortfall/Chelsa"
files_env <- list.files(path_env, pattern = "^climate_Indomalayan.*\\.tif$", full.names = TRUE)
allrasters <- lapply(files_env, raster)
clim_env <- raster::stack(files_env)


# Check the resulting aggregated raster to confirm the new resolution
res(clim_env)
r_IM <- raster(clim_env, 1) #extract 1 raster to check NEW resolution of cells
res_IM <- res(r_IM)  # Save value for later analysis
climate_IndoMalay <- clim_env # create climate stack of our new rasters
plot(climate_IndoMalay[[1]])

# PCA ---------------------------------------------------------------------

library(psych)
library(maps)
library(colorRamps)

# We reduce all the variables to fewer variables using a PCA.
# 1) Here we standardize and prepare the data.
valuesIndoMalay <- as.data.frame(values(climate_IndoMalay)) # get env values from rasters

# 2) Standardization of the variables (Normalization - Mean= 0 and Std= 1)
std <- function(x){(x - mean(x, na.rm = T)) / sd(x, na.rm = T)} 

valuesIndoMalayStd <- apply(valuesIndoMalay, 2, std) # apply to dataframe of climate values
rem <- apply(is.na(valuesIndoMalayStd), 1, any)
PCAIndoMalay <- as.data.frame(valuesIndoMalayStd[!rem, ])

# 3) Run the PCA (2 axis)
matIndoMalay <- matrix(runif(nrow(PCAIndoMalay) * ncol(PCAIndoMalay), 0.00001, 0.00009), 
                       ncol = ncol(PCAIndoMalay)) # add a very small randomness to avoid singularity
PCAIndoMalay2 <- PCAIndoMalay + matIndoMalay
myPCAIndoMalay <- principal(PCAIndoMalay2,
                            nfactors = 2,
                            rotate = "varimax",
                            scores = T)

# 4)
loadings_IM <- as.data.frame(unclass(myPCAIndoMalay$loadings))

loadings_PC1_IM <- loadings_IM %>%
  mutate(abs_RC1 = abs(RC1)) %>%
  arrange(desc(abs_RC1))

loadings_PC2_IM <- loadings_IM %>%
  mutate(abs_RC2 = abs(RC2)) %>%
  arrange(desc(abs_RC2))

# 5) Plot PCA
biplot.psych(myPCAIndoMalay, xlim.s = c(-4, 3), ylim.s = c(-6, 3)) # change limits to plot the whole pca values
rownames(loadings_IM) <- str_extract(rownames(loadings_IM), "bio\\d+")
loadings_IM <- tibble::rownames_to_column(loadings_IM, var = "variables")
scores_IM <- myPCAIndoMalay$scores
scale_factor <- 2.5

biplot_IM <- ggplot() +
    geom_point(data = as.data.frame(scores_IM), aes(x = RC1, y = RC2), color = "black", size = 1) +
    geom_text(data = loadings_IM, aes(x = RC1  * scale_factor * 1.05, y = RC2  * scale_factor * 1.05, label = variables), 
              vjust = -0.5, hjust = 0.5, size = 2, color = "#F0A202", fontface = "bold") +
    geom_segment(data = loadings_IM, aes(x = 0, y = 0, xend = RC1 * scale_factor, yend = RC2 * scale_factor), 
                 arrow = arrow(type = "open", length = unit(0.1, "inches")), color = "#0072B2") +
    scale_x_continuous(breaks = seq(-4, 2, by = 1)) + 
    scale_y_continuous(breaks = seq(-5, 3, by = 1)) +
    xlab(paste("PC1 (39%)")) +
    ylab(paste("PC2 (31%)")) +
    theme_test(base_size = 12) +
    theme(aspect.ratio = 1)

biplot_IM

# PCA of our study area ---------------------------------------------------

# Instead of assuming the absolute values of all Chelsa variables, 
# the next map assumes values of the linear combinations between these variables (PCA scores).
valuesIndoMalay_2 <- valuesIndoMalayStd[, 1:2] # create a vector with the same length as v2 but 2 columns
valuesIndoMalay_2[!rem, ] <- myPCAIndoMalay$scores # insert PCA_scores (2axis-2columns) to the new vector
climate_PCA_IndoMalay <- subset(climate_IndoMalay, 1:2) # extract 2 raster layer as a layerbase to insert our pca values
values(climate_PCA_IndoMalay) <- valuesIndoMalay_2 # insert pca values of v3 into raster
names(climate_PCA_IndoMalay) <- c("PC1", "PC2") # rename
# coords_PCA_IM <- as.data.frame(climate_PCA_IndoMalay, xy = T)
# write.table(coords_PCA_IM, "Processed_Data/coords_PCA_IM.txt") # Save for environmental distance calculation


# Plot study area by each axis from PCA values
world_sf <- ne_countries(scale = "medium", returnclass = "sf") %>%
  dplyr::filter(name != "Antarctica")

par(mfrow = c(1, 2))
plot(climate_PCA_IndoMalay, 1, font = 1, font.lab = 1, add = T)  # PCA axis 1
plot(world_sf, add = TRUE, border = "gray30", col = NA, lwd = 0.8, xlim = c(-12.5, 36), ylim = c(65.85, 134))
plot(climate_PCA_IndoMalay, 2, font = 1, font.lab = 1)  # PCA axis 2
plot(world_sf, add = TRUE, border = "gray30", col = NA, lwd = 0.8, xlim = c(-12.5, 36), ylim = c(65.85, 134))


# Environmental space -----------------------------------------------------

# The next step is to create the environmental space using the two PCA scores.
# Create a Cartesian plan with PC1 and PC2 scores

climate_PCA_values_IndoMalay <- values(climate_PCA_IndoMalay) # get env values from PCA
# Transform the two vars we want into the env space, 
# by taking the min and max scores (PCA score values) for each PCA axis
# and create a raster object. 
xminIM <- min(climate_PCA_values_IndoMalay[, 1], na.rm = TRUE) 
xmaxIM <- max(climate_PCA_values_IndoMalay[, 1], na.rm = TRUE)
yminIM <- min(climate_PCA_values_IndoMalay[, 2], na.rm = TRUE)
ymaxIM <- max(climate_PCA_values_IndoMalay[, 2], na.rm = TRUE)

# This function creates the cartesian plan comprising the min and max PCA scores
env_space_IM <- raster(xmn = xminIM, xmx = xmaxIM, 
                       ymn = yminIM, ymx = ymaxIM, crs = 4326,
                       res = 0.5) # Resolution of the env cell can change

extent(env_space_IM) <- extent(climate_PCA_values_IndoMalay)

ncell(env_space_IM)
values(env_space_IM) <- 0
env_space_area_IM <- env_space_IM # duplicate this object for the next step

### Insert our PCA values into this env_space
# Extract PC1 and PC2, convert it in "classes of values" to plot in the map
env_space_IM_v <- raster::extract(env_space_IM, 
                                  y = na.omit(climate_PCA_values_IndoMalay), 
                                  cellnumbers = TRUE)[, 1]
n_env_space_IM_v <- table(env_space_IM_v) # this function counts the frequency of "climates" (PC scores)

values(env_space_area_IM)[as.numeric(names(n_env_space_IM_v))] <- n_env_space_IM_v
area_values_IM <- values(env_space_area_IM)
area_values_IM[area_values_IM == 0] <- NA
values(env_space_area_IM) <- area_values_IM
pol <- rasterToPolygons(env_space_area_IM, dissolve = FALSE)

# The following map shows the frequency of climates in the study area.
# The axes shows all the PC1 and PC2 scores. 
# Each cell "roughly" represents a climate type.
# The colors indicate the frequency of climate types.

cols <- colorRampPalette(c("#56B4E9", "#33A02C", "#FDE725FF", "#CC79A7", "#D55E00"))

plot(env_space_area_IM, 
     xlab = "PC1", 
     ylab = "PC2",
     col = cols(1000),
     xlim = c(-3.68, 1.83),
     ylim = c(-4.57, 2.49),
     cex.axis = 0.8,
     cex.lab = 0.8,
     axis.args = list(cex.axis = 0.8),
     font = 2, font.lab = 2
)
plot(pol, add = TRUE)

# writeRaster(env_space_area_IM, "Processed_Data/env_space_IM.tif", format = "GTiff", overwrite = TRUE)

# Well-sampled cells ------------------------------------------------------

# Import well-sampled cells with a slope threshold of 0.1 from the inventory completeness analysis
# well_sampled_2 <- read.csv("Data/well_sampled_2.csv")

# Select only well-sampled cells within the boundaries of Indomalayan realm
well_sampled_IM_2 <- well_sampled_2 %>% 
  dplyr::filter(Realm == "IM")

#_______________________________________________________________________________
# NOW we calculate the environmental space for our Well-sampled cells
# Using estimators output from 'InventoryCompleteness' script:
well_sampled_IM_2 <- well_sampled_IM_2[, c('Longitude','Latitude')]

# Plot the well-sampled occurrences on the map
plot(r_IM, main = "", col = "gray", legend = FALSE)
points(well_sampled_IM_2, col = rgb(1, 0, 0, .5), pch = 20, cex = 2)

# Extract climate values of Well-sampled cells
values_WS_IM <- raster::extract(climate_IndoMalay, 
                                well_sampled_IM_2,
                                cellnumbers = TRUE)[, 1]
# write.table(values_WS_IM, "Processed_Data/values_WS_IM.txt") # Save for environmental distance calculation

coords_WS_IM <- climate_PCA_values_IndoMalay[values_WS_IM, ]

# write.table(coords_WS_IM, "Processed_Data/coords_WS_IM.txt") # Save for environmental distance calculation
cell_WS_IM <- raster::extract(env_space_IM, 
                              coords_WS_IM,
                              cellnumbers = TRUE)[, 1]
n_WS_IM <- table(cell_WS_IM)
env_space_WS_IM <- env_space_IM
values(env_space_WS_IM)[as.numeric(names(n_WS_IM))] <- n_WS_IM

# The following map shows the frequency of climates in the well-sampled cells 
# in the study area (legend of this figure is continuous BUT values are integer numbers, 
# so manually adapt the legend to each case)

plot(env_space_WS_IM, 
     xlab = "PC1", 
     ylab = "PC2",
     col = c("transparent", cols(1000)),
     xlim = c(-3.68, 1.83),
     ylim = c(-4.57, 2.49),
     cex.axis = 0.8,
     cex.lab = 0.8,
     axis.args = list(cex.axis = 0.8),
     font = 2, font.lab = 2)
plot(pol, add = TRUE)

# writeRaster(env_space_WS_IM, "Processed_Data/env_space_WS_IM.tif", format = "GTiff", overwrite = TRUE)

# Schoener's D ------------------------------------------------------------

# Schoener's D: quantifies the overlap between the location of well-sampled 
# sites and cells with most frequent conditions 
# The Schoener's D index varies from zero (total lack of congruence) to one 
# (total congruence) 

# Transform the abundance of each cell into probabilities.
area_values_IM <- area_values_IM/sum(area_values_IM, na.rm = TRUE) # Relative frequency of climate type for all the study area
WS_values_IM <- values(env_space_WS_IM)
WS_values_IM <- WS_values_IM/sum(WS_values_IM, na.rm = TRUE) # Relative frequency of climate type for well-sampled cells

# Calculate the climate overlap using Schoener's D.
# D values close to 1 indicate that the location of well-sampled 
# sites coincide with climate conditions 
# that are frequently found in the study area
SchoenersD <- function(x, y) {
  sub_values <- abs(x - y)
  D <- 1 - (sum(sub_values, na.rm = TRUE) / 2)
  return(D)
}

D_IM <- SchoenersD(area_values_IM, WS_values_IM)
print(paste("Climate overlap between well-sampled cells and the study area, given by the observed Schoener's D equals = ", 
            round(D_IM, 3), "%"))

# We create a null model to test if D values is different from a 
# random distribution of D values calculated from randomly sampling occurrence 
# records.
set.seed(0)
replications <- 1000 # Choose the number of replications
D_rnd <- numeric(replications)
for (i in 1:replications) {
  rnd <- sample(env_space_IM_v, length(well_sampled_IM_2), replace = TRUE)
  n_rnd <- table(rnd)
  env_space_rnd <- env_space_IM
  values(env_space_rnd)[as.numeric(names(n_rnd))] <- n_rnd
  rnd_values <- values(env_space_rnd)
  rnd_values <- rnd_values/sum(rnd_values, na.rm = TRUE)
  D_rnd[i] <- SchoenersD(area_values_IM, rnd_values)
}

# p-value from previous analysis
# If p <0.05, it means that the location of well-sampled sites does not 
# coincide with areas with climate conditions frequently found in your study area
p_IM <- (sum(D_IM > D_rnd) + 1) / (length(D_rnd) + 1) # Unicaudal test
print(paste("p value equals = ", round(p_IM, 3)))

# Plot distribution of model values
hist(D_rnd, 10, # Write number of bins
     xlim = c(0, 1), 
     main = (""),
     xlab = "D", 
     col = rgb(.5, .5, .5), 
     border = FALSE,
     font = 1, font.lab = 1)
abline(v = D_IM, col = "red", lty = 2) # Red line show the observed Schoeners' D value

############ Kruskal-Wallis test ##################################
# The following map shows the distribution of well-sampled cells (red bars) 
# vs the study area cells (gray bars)
### For a) PC1
hist(myPCAIndoMalay$scores[, 1], 
     breaks = ncol(env_space_IM),
     freq = F,
     col = "grey", 
     border = FALSE, 
     xlim= c(-4, 3), 
     ylim = c(0, 0.5), 
     main = "",
     xlab = "PC1", 
     font = 2, font.lab = 2, 
     cex.lab = 0.8, cex.axis = 0.8)
lines(density(na.omit(myPCAIndoMalay$scores[, 1])), col = "black", lwd = 2)

hist(coords_WS_IM[, 1],
     add = TRUE,
     col = adjustcolor("#0072B2", alpha.f = 0.5),
     freq = F,
     border = FALSE)
lines(density(na.omit(coords_WS_IM[, 1])), col = "#0072B2", lwd = 2)

# X axis is a probability density
# Kruskal-Wallis verifies whether 1) the distribution of well-sampled sites 
# is an unbiased subset of the entire climate conditions of the Indo-Malay. 
# If this is so, p > 0.05 
x_1_IM <- c(myPCAIndoMalay$scores[, 1], coords_WS_IM[, 1])
g_1_IM <- as.factor(c(rep("area", length(myPCAIndoMalay$scores[, 1])),
                      rep("WS", length(coords_WS_IM[, 1]))))
kruskal.test(x_1_IM ~ g_1_IM)
# Kolmogorov smirnov test###
ks.test(myPCAIndoMalay$scores[, 1], coords_WS_IM[, 1])

### For b) PC2
hist(myPCAIndoMalay$scores[, 2], 
     breaks = ncol(env_space_IM),
     freq = F,
     col = "grey", 
     border = FALSE, 
     xlim = c(-6, 4), 
     ylim = c(0, 0.7), 
     main = "",
     xlab = "PC2", 
     font = 2, font.lab = 2, 
     cex.lab = 0.8, cex.axis = 0.8)
lines(density(na.omit(myPCAIndoMalay$scores[, 2])), col = "black", lwd = 2)

hist(coords_WS_IM[, 2],
     breaks = ncol(env_space_IM),
     add = TRUE,
     col = adjustcolor("#0072B2", alpha.f = 0.5),
     freq = F,
     border = FALSE)
lines(density(na.omit(coords_WS_IM[, 2])), col = "#0072B2", lwd = 2)

x_2_IM <- c(myPCAIndoMalay$scores[, 2], coords_WS_IM[, 2])
g_2_IM <- as.factor(c(rep("area", length(myPCAIndoMalay$scores[, 2])),
                      rep("WS", length(coords_WS_IM[, 2]))))
kruskal.test(x_2_IM ~ g_2_IM)
# Kolmogorov smirnov test###
ks.test(myPCAIndoMalay$scores[, 2], coords_WS_IM[, 2])

# Rarity ------------------------------------------------------------------

# We can check how many environmental space has been sampled 
# and how does these cells look like.

surface_IM <- WS_values_IM
surface_IM[is.na(area_values_IM) | area_values_IM == 0] <- NA
sampled_IM <- sum(surface_IM > 0, na.rm = TRUE)/ sum(surface_IM >= 0, na.rm = TRUE) * 100
percen_IM <- round(sampled_IM, 2)
prin_IM <- paste0(percen_IM, "% of our study area climate types covered by well-sampled cells")
print(prin_IM)

# Is this env. space sampled corresponding to rare climates?
# First, we make values vary from 0 to 1 according to their rarity:
# Values close to 0, are very common, values close to 1 very rare
mini_IM <- min(area_values_IM, na.rm = TRUE) # less frequent value

# rarity index, also called Min-Max scalling
area_values01_IM <- abs(1 - (area_values_IM - mini_IM) / (max(area_values_IM, na.rm = TRUE) - mini_IM)) 
# see http://rasbt.github.io/mlxtend/user_guide/preprocessing/minmax_scaling/
hist(area_values01_IM, 
     breaks = ncol(env_space_IM), 
     freq = FALSE, col = "white",
     main = "", 
     xlab = "Climate rarity", border = FALSE, 
     ylim = c(0, 5),
     font = 2, font.lab = 2,
     cex.lab = 0.8, cex.axis = 0.8)
lines(density(na.omit(area_values01_IM)), col = "black", lwd = 2)
lines(density(na.omit(area_values01_IM[surface_IM > 0])), col = "#0072B2", lwd = 2)

# Kruskal-Wallis test to see if the distribution of rarities for the sampled 
# occurrences differs from the distribution observed in the entire area
x_IM <- c(area_values01_IM, area_values01_IM[surface_IM > 0])
g_IM <- as.factor(c(rep("area", length(area_values01_IM)),
                    rep("ws", length(area_values01_IM[surface_IM > 0]))))
kruskal.test(x_IM ~ g_IM)
# Kolmogorov smirnov test
area_values03_IM <- (na.omit(area_values01_IM))
area_values04_IM <- (na.omit(area_values01_IM[surface_IM > 0]))
ks.test(area_values03_IM, area_values04_IM)

# This result indicates that the under-sampled area is composed 
# mainly by rare climates. You can inspect the environmental space figures 
# to check which areas were not sampled.

# Finally MAP the climatic rarity
rarity_env_IM <- env_space_IM
values(rarity_env_IM) <- area_values01_IM
rarity_Percell_IM <- raster::extract(rarity_env_IM, climate_PCA_values_IndoMalay)
rarity_map_IM <- r_IM
values(rarity_map_IM) <- rarity_Percell_IM

cols <- colorRampPalette(c('#e0f3db','#a8ddb5','#43a2ca'))
x11()
plot(rarity_map_IM, 
     col = cols(10), 
     font = 2, font.lab = 2, 
     ylim = c(-12.19911, 35.82528), 
     xlim = c(65.85641, 131.0855)
)

points(well_sampled_IM_2, col = rgb(1, 0, 0, .5), pch = 19, cex = 1)

# Convert raster to data frame
rarity_map_df_IM <- as.data.frame(as(rarity_map_IM, "SpatialPixelsDataFrame"))
colnames(rarity_map_df_IM) <- c("value", "x", "y")

# write.csv(rarity_map_df_IM, "Processed_Data/Rarity_IM.csv") # Save to create rarity map



