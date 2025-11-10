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

# Import the 19 CHELSA bioclimatic variables (1° resolution) cropped to the Nearctic realm

#Raster

path_env <- "D:/Scrutinizing-the-Wallacean-shortfall/Chelsa"
files_env <- list.files(path_env, pattern = "^climate_Nearctic.*\\.tif$", full.names = TRUE)
allrasters <- lapply(files_env, raster)
clim_env <- raster::stack(files_env)

# Check the resulting aggregated raster to confirm the new resolution
res(clim_env)
r_NA <- raster(clim_env, 1) #extract 1 raster to check NEW resolution of cells
res_NA <- res(r_NA)  # Save value for later analysis
climate_Nearctic <- clim_env # create climate stack of our new rasters
plot(climate_Nearctic[[2]])


# PCA ---------------------------------------------------------------------

library(psych)
library(maps)
library(colorRamps)

# We reduce all the variables to fewer variables using a PCA.
# 1) Here we standardize and prepare the data.
valuesNearctic <- as.data.frame(values(climate_Nearctic)) # get env values from rasters

# 2) Standardization of the variables (Normalization - Mean = 0 and Std = 1)
std <- function(x){(x - mean(x, na.rm = T)) / sd(x, na.rm = T)} 
valuesNearcticStd <- apply(valuesNearctic, 2, std) # apply to dataframe of climate values

rem <- apply(is.na(valuesNearcticStd), 1, any)
PCANearctic <- as.data.frame(valuesNearcticStd[!rem, ])

# 3) Run the PCA (2 axis)
matNearctic <- matrix(runif(nrow(PCANearctic) * ncol(PCANearctic), 0.00001, 0.00009), 
                      ncol = ncol(PCANearctic)) # add a very small randomness to avoid singularity
PCANearctic2 <- PCANearctic + matNearctic
myPCANearctic <- principal(PCANearctic2,
                           nfactors = 2,
                           rotate = "varimax",
                           scores = T)

# 4)
loadings_NA <- as.data.frame(unclass(myPCANearctic$loadings))

loadings_PC1_NA <- loadings_NA %>%
  mutate(abs_RC1 = abs(RC1)) %>%
  arrange(desc(abs_RC1))

loadings_PC2_NA <- loadings_NA %>%
  mutate(abs_RC2 = abs(RC2)) %>%
  arrange(desc(abs_RC2))


# 5) Plot PCA
biplot.psych(myPCANearctic, xlim.s = c(-4, 4), ylim.s = c(-2, 3))
rownames(loadings_NA) <- str_extract(rownames(loadings_NA), "bio\\d+")
loadings_NA <- tibble::rownames_to_column(loadings_NA, var = "variables")
scores_NA <- myPCANearctic$scores
scale_factor <- 2.5

biplot_NA <- ggplot() +
    geom_point(data = as.data.frame(scores_NA), aes(x = RC1, y = RC2), color = "black", size = 1) +
    geom_text(data = loadings_NA, aes(x = RC1  * scale_factor * 1.05, y = RC2  * scale_factor * 1.05, label = variables), 
              vjust = -0.5, hjust = 0.5, size = 2, color = "#F0A202", fontface = "bold") +
    geom_segment(data = loadings_NA, aes(x = 0, y = 0, xend = RC1 * scale_factor, yend = RC2 * scale_factor), 
                 arrow = arrow(type = "open", length = unit(0.1, "inches")), color = "#0072B2") +
    scale_x_continuous(breaks = seq(-2, 8, by = 2)) + 
    scale_y_continuous(breaks = seq(-2, 3, by = 2)) +
    xlab(paste("PC1 (38%)")) +
    ylab(paste("PC2 (36%)")) +
    theme_test(base_size = 12) +
    theme(aspect.ratio = 1)

biplot_NA

# PCA of our study area ---------------------------------------------------

# Instead of assuming the absolute values of all Chelsa variables, 
# the next map assumes values of the linear combinations between these variables (PCA scores).
valuesNearctic_2 <- valuesNearcticStd[, 1:2] # create a vector with the same length as v2 but 2 columns
valuesNearctic_2[!rem, ] <- myPCANearctic$scores 
climate_PCA_Nearctic <- subset(climate_Nearctic, 1:2)
values(climate_PCA_Nearctic) <- valuesNearctic_2 # insert pca values of v3 into raster
names(climate_PCA_Nearctic) <- c("PC1", "PC2") # rename
# coords_PCA_NA <- as.data.frame(climate_PCA_Nearctic, xy = T)
# write.table(coords_PCA_NA, "Processed_Data/coords_PCA_NA.txt") # Save for environmental distance calculation

# Plot study area by each axis from PCA values
world_sf <- ne_countries(scale = "medium", returnclass = "sf") %>%
  dplyr::filter(name != "Antarctica")

par(mfrow = c(1, 2))
plot(climate_PCA_Nearctic, 1, font = 1, font.lab = 1, add = T, xlim = c(-172, -18), ylim = c(20, 85))  # PCA axis 1
plot(world_sf, add = TRUE, border = "gray30", col = NA, lwd = 0.8)
plot(climate_PCA_Nearctic, 2, font = 1, font.lab = 1, xlim = c(-172, -18), ylim = c(20, 85))  # PCA axis 2
plot(world_sf, add = TRUE, border = "gray30", col = NA, lwd = 0.8)


# Environmental space -----------------------------------------------------

# The next step is to create the environmental space using the two PCA scores.
# Create a Cartesian plan with PC1 and PC2 scores

climate_PCA_values_Nearctic <- values(climate_PCA_Nearctic) # get env values from PCA

# Transform the two vars we want into the env space, 
# by taking the min and max scores (PCA score values) for each PCA axis
# and create a raster object. 
xminNA <- min(climate_PCA_values_Nearctic[, 1], na.rm = TRUE) 
xmaxNA <- max(climate_PCA_values_Nearctic[, 1], na.rm = TRUE)
yminNA <- min(climate_PCA_values_Nearctic[, 2], na.rm = TRUE)
ymaxNA <- max(climate_PCA_values_Nearctic[, 2], na.rm = TRUE)

# This function creates the cartesian plan comprising the min and max PCA scores

env_space_NA <- raster(xmn = xminNA, xmx = xmaxNA, 
                       ymn = yminNA, ymx = ymaxNA, crs = 4326,
                       res = 0.5) # Resolution of the env cell can change
extent(env_space_NA) <- extent(climate_PCA_values_Nearctic)
ncell(env_space_NA)
values(env_space_NA) <- 0
env_space_area_NA <- env_space_NA # duplicate this object for the next step

### Insert our PCA values into this env_space
# Extract PC1 and PC2, convert it in "classes of values" to plot in the map
env_space_NA_v <- raster::extract(env_space_NA, 
                                  y = na.omit(climate_PCA_values_Nearctic), 
                                  cellnumbers = TRUE)[, 1]
table(env_space_NA_v)
n_env_space_NA_v <- table(env_space_NA_v) # this function counts the frequency of "climates" (PC scores)


values(env_space_area_NA)[as.numeric(names(n_env_space_NA_v))] <- n_env_space_NA_v
area_values_NA <- values(env_space_area_NA)
area_values_NA[area_values_NA == 0] <- NA
values(env_space_area_NA) <- area_values_NA
pol <- rasterToPolygons(env_space_area_NA, dissolve = FALSE)

# The following map shows the frequency of climates in the study area.
# The axes shows all the PC1 and PC2 scores. 
# Each cell "roughly" represents a climate type.
# The colors indicate the frequency of climate types.

cols <- colorRampPalette(c("#56B4E9", "#33A02C", "#FDE725FF", "#CC79A7", "#D55E00"))

plot(env_space_area_NA, 
     xlab = "PC1", 
     ylab = "PC2",
     col = cols(1000),
     xlim = c(-1.52, 8),
     ylim = c(-1.92, 2.44),
     cex.axis = 0.8,
     cex.lab = 0.8,
     axis.args = list(cex.axis = 0.8),
     font = 2, font.lab = 2
)
plot(pol, add = TRUE)

# writeRaster(env_space_area_NA, "Processed_Data/env_space_NA.tif", format = "GTiff", overwrite = TRUE)


# Well-sampled cells ------------------------------------------------------

# Import well-sampled cells with a slope threshold of 0.05 from the inventory completeness analysis
# well_sampled_1 <- read.csv("Data/well_sampled_1.csv")
well_sampled_1 <- well_sampled_1 %>%
  mutate(Realm = ifelse(is.na(Realm), "NA", Realm))

# Select only well-sampled cells within the boundaries of Nearctic realm
well_sampled_NA_1 <- well_sampled_1 %>% 
  dplyr::filter(Realm == "NA")

#_______________________________________________________________________________
# Now we calculate the environmental space for our well-sampled cells
# Using estimators output from 'InventoryCompleteness' script:
well_sampled_NA_1 <- well_sampled_NA_1[, c('Longitude','Latitude')]

# Plot the well-sampled occurrences on the map
raster::plot(r_NA, main = "", col = "gray", legend = FALSE)
points(well_sampled_NA_1, col = rgb(1, 0, 0, .5), pch = 20, cex = 2)

# Extract climate values of Well-sampled cells
values_WS_NA <- raster::extract(climate_Nearctic, 
                                well_sampled_NA_1,
                                cellnumbers = TRUE)[, 1]

# write.table(values_WS_NA, "Processed_Data/values_WS_NA.txt") # Save for environmental distance calculation

coords_WS_NA <- climate_PCA_values_Nearctic[values_WS_NA, ]
# write.table(coords_WS_NA, "Processed_Data/coords_WS_NA.txt") # Save for environmental distance calculation

cell_WS_NA <- raster::extract(env_space_NA, 
                              coords_WS_NA,
                              cellnumbers = TRUE)[, 1]

n_WS_NA <- table(cell_WS_NA)
env_space_WS_NA <- env_space_NA
values(env_space_WS_NA)[as.numeric(names(n_WS_NA))] <- n_WS_NA
ncell(cell_WS_NA)

# The following map shows the frequency of climates in the well-sampled cells 
# in the study area (legend of this figure is continuous BUT values are integer numbers, 
# so manually adapt the legend to each case)

plot(env_space_WS_NA, 
     xlab = "PC1", 
     ylab = "PC2",
     col = c("transparent", cols(1000)),
     xlim = c(-1.52, 8),
     ylim = c(-1.92, 2.44),
     cex.axis = 0.8,
     cex.lab = 0.8,
     axis.args = list(cex.axis = 0.8), 
     font = 2, font.lab = 2)
plot(pol, add = TRUE)

# writeRaster(env_space_WS_NA, "Processed_Data/env_space_WS_NA.tif", format = "GTiff", overwrite = TRUE)
# Schoener's D ------------------------------------------------------------

# Schoener's D: quantifies the overlap between the location of well-sampled 
# sites and cells with most frequent conditions 
# The Schoener's D index varies from zero (total lack of congruence) to one 
# (total congruence) 

# Transform the abundance of each cell into probabilities.
area_values_NA <- area_values_NA/sum(area_values_NA, na.rm = TRUE) # Relative frequency of climate type for all the study area
WS_values_NA <- values(env_space_WS_NA)
WS_values_NA <- WS_values_NA/sum(WS_values_NA, na.rm = TRUE) # Relative frequency of climate type for well-sampled cells

# Calculate the climate overlap using Schoener's D.
# D values close to 1 indicate that the location of well-sampled 
# sites coincide with climate conditions 
# that are frequently found in the study area
SchoenersD <- function(x, y) {
  sub_values <- abs(x - y)
  D <- 1 - (sum(sub_values, na.rm = TRUE) / 2)
  return(D)
}

D_NA <- SchoenersD(area_values_NA, WS_values_NA)
print(paste("Climate overlap between well-sampled cells and the study area, given by the observed Schoener's D equals = ", 
            round(D_NA,3), "%"))

# We create a null model to test if D values is different from a 
# random distribution of D values calculated from randomly sampling occurrence 
# records.
set.seed(0)
replications <- 1000 # Choose the number of replications
D_rnd <- numeric(replications)
for (i in 1:replications) {
  rnd <- sample(env_space_NA_v, length(well_sampled_NA_1), replace = TRUE)
  n_rnd <- table(rnd)
  env_space_rnd <- env_space_NA
  values(env_space_rnd)[as.numeric(names(n_rnd))] <- n_rnd
  rnd_values <- values(env_space_rnd)
  rnd_values <- rnd_values/sum(rnd_values, na.rm = TRUE)
  D_rnd[i] <- SchoenersD(area_values_NA, rnd_values)
}

# p-value from previous analysis
# If p < 0.05, it means that the location of well-sampled sites does not 
# coincide with areas with climate conditions frequently found in your study area
p_NA <- (sum(D_NA > D_rnd) + 1) / (length(D_rnd) + 1) # Unicaudal test
print(paste("p value equals = ", round(p_NA, 3)))

# Plot distribution of model values
hist(D_rnd, 10, # Write number of bins
     xlim = c(0, 1), 
     main = (""),
     xlab = "D", 
     col = rgb(.5, .5, .5), 
     border = FALSE,
     font = 1, font.lab = 1)
abline(v = D_NA, col = "red", lty = 2) # Red line show the observed Schoeners' D value

############ Kruskal-Wallis test ##################################
# The following map shows the distribution of well-sampled cells (red bars) 
# vs the study area cells (gray bars)
### For a) PC1

hist(myPCANearctic$scores[, 1], 
     breaks = ncol(env_space_NA),
     freq = F,
     col = "grey", 
     border = FALSE, 
     xlim= c(-2, 8), 
     ylim = c(0, 0.65), 
     main = "",
     xlab = "PC1", 
     font = 2, font.lab = 2, 
     cex.lab = 0.8, cex.axis = 0.8)
lines(density(na.omit(myPCANearctic$scores[, 1])), col = "black", lwd = 2)

hist(coords_WS_NA[, 1],
     add = TRUE,
     col = adjustcolor("#0072B2", alpha.f = 0.5),
     freq = F,
     border = FALSE)
lines(density(na.omit(coords_WS_NA[, 1])), col = "#0072B2", lwd = 2)

# X axis is a probability density
# Kruskal-Wallis verifies whether 1) the distribution of well-sampled sites 
# is an unbiased subset of the entire climate conditions of the World. 
# If this is so, p > 0.05 
x_1_NA <- c(myPCANearctic$scores[, 1], coords_WS_NA[, 1])
g_1_NA <- as.factor(c(rep("area", length(myPCANearctic$scores[, 1])),
                      rep("WS", length(coords_WS_NA[, 1]))))
kruskal.test(x_1_NA ~ g_1_NA)
# Kolmogorov smirnov test###
ks.test(myPCANearctic$scores[, 1], coords_WS_NA[, 1])

### For b) PC2
hist(myPCANearctic$scores[, 2], 
     breaks = ncol(env_space_NA),
     freq = F,
     col = "grey", 
     border = FALSE, 
     xlim = c(-2.8, 3), 
     ylim = c(0, 0.65), 
     main = "",
     xlab = "PC2", 
     font = 2, font.lab = 2, 
     cex.lab = 0.8, cex.axis = 0.8)
lines(density(na.omit(myPCANearctic$scores[, 2])), col = "black", lwd = 2)

hist(coords_WS_NA[, 2],
     breaks = ncol(env_space_NA), 
     add = TRUE,
     col = adjustcolor("#0072B2", alpha.f = 0.5),
     freq = F,
     border = FALSE)
lines(density(na.omit(coords_WS_NA[, 2])), col = "#0072B2", lwd = 2)

x_2_NA <- c(myPCANearctic$scores[, 2], coords_WS_NA[, 2])
g_2_NA <- as.factor(c(rep("area", length(myPCANearctic$scores[, 2])),
                      rep("WS", length(coords_WS_NA[, 2]))))
kruskal.test(x_2_NA ~ g_2_NA)

# Kolmogorov smirnov test###
ks.test(myPCANearctic$scores[, 2], coords_WS_NA[, 2])

# Rarity ------------------------------------------------------------------

# We can check how many environmental space has been sampled 
# and how does these cells look like.
surface_NA <- WS_values_NA
surface_NA[is.na(area_values_NA) | area_values_NA == 0] <- NA
sampled_NA <- sum(surface_NA > 0, na.rm = TRUE)/ sum(surface_NA >= 0, na.rm = TRUE) * 100
percen_NA <- round(sampled_NA, 2)
prin_NA <- paste0(percen_NA, "% of our study area climate types covered by well-sampled cells")
print(prin_NA)

# Is this env. space sampled corresponding to rare climates?
# First, we make values vary from 0 to 1 according to their rarity:
# Values close to 0, are very common, values close to 1 very rare
mini_NA <- min(area_values_NA, na.rm = TRUE) # less frequent value
# rarity index, also called Min-Max scalling
area_values01NA <- abs(1 - (area_values_NA - mini_NA) / (max(area_values_NA, na.rm = TRUE) - mini_NA)) 
# see http://rasbt.github.io/mlxtend/user_guide/preprocessing/minmax_scaling/
hist(area_values01NA, 
     breaks = ncol(env_space_NA), 
     freq = FALSE, col = "white",
     main = "", 
     xlab = "Climate rarity", border = FALSE, 
     ylim = c(0, 5),
     font = 2, font.lab = 2,
     cex.lab = 0.8, cex.axis = 0.8)
lines(density(na.omit(area_values01NA)), col = "black", lwd = 2)
lines(density(na.omit(area_values01NA[surface_NA > 0])), col = "#0072B2", lwd = 2)

# Kruskal-Wallis test to see if the distribution of rarities for the sampled 
# occurrences differs from the distribution observed in the entire area
x_NA <- c(area_values01NA, area_values01NA[surface_NA > 0])
g_NA <- as.factor(c(rep("area", length(area_values01NA)),
                    rep("ws", length(area_values01NA[surface_NA > 0]))))
kruskal.test(x_NA ~ g_NA)
# Kolmogorov smirnov test
area_values03_NA <- (na.omit(area_values01NA))
area_values04_NA <- (na.omit(area_values01NA[surface_NA > 0]))
ks.test(area_values03_NA, area_values04_NA)

# This result indicates that the under-sampled area is composed 
# mainly by rare climates. You can inspect the environmental space figures 
# to check which areas were not sampled.

# Finally MAP the climatic rarity
rarity_env_NA <- env_space_NA
values(rarity_env_NA) <- area_values01NA
rarity_Percell_NA <- raster::extract(rarity_env_NA, climate_PCA_values_Nearctic)
rarity_map_NA <- r_NA
values(rarity_map_NA) <- rarity_Percell_NA

cols <- colorRampPalette(c('#e0f3db','#a8ddb5','#43a2ca'))
x11()
plot(rarity_map_NA,
     col = cols(10),
     font = 2, font.lab = 2,
     ylim = c(19.10538, 83.62313),
     xlim = c(-179.142, 179.7775)
)

points(well_sampled_NA_1, col = rgb(1, 0, 0, .5), pch = 19, cex = 1)

# Convert raster to data frame 
rarity_map_df_NA <- as.data.frame(as(rarity_map_NA, "SpatialPixelsDataFrame"))
colnames(rarity_map_df_NA) <- c("value", "x", "y")

# write.csv(rarity_map_df_NA, "Processed_Data/Rarity_NA.csv") # Save to create rarity map

