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


# Import the 19 CHELSA bioclimatic variables (1° resolution) cropped to the Eastern Palearctic realm

#Raster
path_env <- "D:/Scrutinizing-the-Wallacean-shortfall/Chelsa"
files_env <- list.files(path_env, pattern = "^climate_E_Palearctic.*\\.tif$", full.names = TRUE)
allrasters <- lapply(files_env, raster)
clim_env <- raster::stack(files_env)

# Check the resulting aggregated raster to confirm the new resolution
res(clim_env)
r_EPA <- raster(clim_env, 1) #extract 1 raster to check NEW resolution of cells
res_EPA <- res(r_EPA)  # Save value for later analysis
climate_E_Palearctic <- clim_env # create climate stack of our new rasters
plot(climate_E_Palearctic[[1]])

# PCA ---------------------------------------------------------------------

library(psych)
library(maps)
library(colorRamps)

# We reduce all the variables to fewer variables using a PCA.
# 1) Here we standardize and prepare the data.
values_E_Palearctic <- as.data.frame(values(climate_E_Palearctic)) # get env values from rasters

# 2) Standardization of the variables (Normalization - Mean = 0 and Std = 1)
std <- function(x){(x - mean(x, na.rm = T)) / sd(x, na.rm = T)} 
valuesEPalearcticStd <- apply(values_E_Palearctic, 2, std) # apply to dataframe of climate values

rem <- apply(is.na(valuesEPalearcticStd), 1, any)
PCAEPalearctic <- as.data.frame(valuesEPalearcticStd[!rem, ])

# 3) Run the PCA (2 axis)
matEPalearctic <- matrix(runif(nrow(PCAEPalearctic) * ncol(PCAEPalearctic), 0.00001, 0.00009), 
                         ncol = ncol(PCAEPalearctic)) # add a very small randomness to avoid singularity
PCAEPalearctic2 <- PCAEPalearctic + matEPalearctic
myPCAEPalearctic <- principal(PCAEPalearctic2,
                              nfactors = 2,
                              rotate = "varimax",
                              scores = T)

# 4)
loadings_EP <- as.data.frame(unclass(myPCAEPalearctic$loadings))

loadings_PC1_EP <- loadings_EP %>%
  mutate(abs_RC1 = abs(RC1)) %>%
  arrange(desc(abs_RC1))

loadings_PC2_EP <- loadings_EP %>%
  mutate(abs_RC2 = abs(RC2)) %>%
  arrange(desc(abs_RC2))

# 5) Plot PCA
biplot.psych(myPCAEPalearctic, xlim.s = c(-2, 3.2), ylim.s = c(-2, 8)) 
rownames(loadings_EP) <- str_extract(rownames(loadings_EP), "bio\\d+")
loadings_EP <- tibble::rownames_to_column(loadings_EP, var = "variables")
scores_EP <- myPCAEPalearctic$scores
scale_factor <- 2.5

biplot_EP <- ggplot() +
    geom_point(data = as.data.frame(scores_EP), aes(x = RC1, y = RC2), color = "black", size = 1) +
    geom_text(data = loadings_EP, aes(x = RC1  * scale_factor * 1.05, y = RC2  * scale_factor * 1.05, label = variables), 
              vjust = -0.5, hjust = 0.5, size = 2, color = "#F0A202", fontface = "bold") +
    geom_segment(data = loadings_EP, aes(x = 0, y = 0, xend = RC1 * scale_factor, yend = RC2 * scale_factor), 
                 arrow = arrow(type = "open", length = unit(0.1, "inches")), color = "#0072B2") +
    scale_x_continuous(breaks = seq(-2, 3, by = 2)) + 
    scale_y_continuous(breaks = seq(-2, 7, by = 2)) +
    xlab(paste("PC1 (34%)")) +
    ylab(paste("PC2 (32%)")) +
    theme_test(base_size = 12) +
    theme(aspect.ratio = 1)

biplot_EP

# PCA of our study area ---------------------------------------------------

# Instead of assuming the absolute values of all Chelsa variables, 
# the next map assumes values of the linear combinations between these variables (PCA scores).
valuesEPalearctic_2 <- valuesEPalearcticStd[, 1:2] # create a vector with the same length as v2 but 2 columns
valuesEPalearctic_2[!rem, ] <- myPCAEPalearctic$scores 
climate_PCA_EPalearctic <- subset(climate_E_Palearctic, 1:2)
values(climate_PCA_EPalearctic) <- valuesEPalearctic_2 # insert pca values of v3 into raster
names(climate_PCA_EPalearctic) <- c("PC1", "PC2") # rename
# coords_PCA_EP <- as.data.frame(climate_PCA_EPalearctic, xy = T)
# write.table(coords_PCA_EP, "Processed_Data/coords_PCA_EP.txt") # Save for environmental distance calculation


# Plot study area by each axis from PCA values
world_sf <- ne_countries(scale = "medium", returnclass = "sf") %>%
  dplyr::filter(name != "Antarctica")

par(mfrow = c(1, 2))
plot(climate_PCA_EPalearctic, 1, font = 1, font.lab = 1, add = T, xlim = c(45, 174), ylim = c(24.5, 80))  # PCA axis 1
plot(world_sf, add = TRUE, border = "gray30", col = NA, lwd = 0.8, xlim = c(45, 174), ylim = c(24.5, 80))
plot(climate_PCA_EPalearctic, 2, font = 1, font.lab = 1, xlim = c(45, 174), ylim = c(24.5, 80))  # PCA axis 2
plot(world_sf, add = TRUE, border = "gray30", col = NA, lwd = 0.8, xlim = c(45, 174), ylim = c(24.5, 80))

# Environmental space -----------------------------------------------------

# The next step is to create the environmental space using the two PCA scores.
# Create a Cartesian plan with PC1 and PC2 scores

climate_PCA_values_EPalearctic <- values(climate_PCA_EPalearctic) # get env values from PCA

# Transform the two vars we want into the env space, 
# by taking the min and max scores (PCA score values) for each PCA axis
# and create a raster object. 
xminEPA <- min(climate_PCA_values_EPalearctic[, 1], na.rm = TRUE) 
xmaxEPA <- max(climate_PCA_values_EPalearctic[, 1], na.rm = TRUE)
yminEPA <- min(climate_PCA_values_EPalearctic[, 2], na.rm = TRUE)
ymaxEPA <- max(climate_PCA_values_EPalearctic[, 2], na.rm = TRUE)

# This function creates the cartesian plan comprising the min and max PCA scores

env_space_EPA <- raster(xmn = xminEPA, xmx = xmaxEPA, 
                        ymn = yminEPA, ymx = ymaxEPA, crs = 4326,
                        res = 0.5) # Resolution of the env cell can change

extent(env_space_EPA) <- extent(climate_PCA_values_EPalearctic)

ncell(env_space_EPA)
values(env_space_EPA) <- 0
env_space_area_EPA <- env_space_EPA # duplicate this object for the next step

### Insert our PCA values into this env_space
# Extract PC1 and PC2, convert it in "classes of values" to plot in the map
env_space_EPA_v <- raster::extract(env_space_EPA, 
                                   y = na.omit(climate_PCA_values_EPalearctic), 
                                   cellnumbers = TRUE)[, 1]
table(env_space_EPA_v)
n_env_space_EPA_v <- table(env_space_EPA_v) # this function counts the frequency of "climates" (PC scores)

values(env_space_area_EPA)[as.numeric(names(n_env_space_EPA_v))] <- n_env_space_EPA_v
area_values_EPA <- values(env_space_area_EPA)
area_values_EPA[area_values_EPA == 0] <- NA
values(env_space_area_EPA) <- area_values_EPA
pol <- rasterToPolygons(env_space_area_EPA, dissolve = FALSE)

# The following map shows the frequency of climates in the study area.
# The axes shows all the PC1 and PC2 scores. 
# Each cell "roughly" represents a climate type.
# The colors indicate the frequency of climate types.

cols <- colorRampPalette(c("#56B4E9", "#33A02C", "#FDE725FF", "#CC79A7", "#D55E00"))

plot(env_space_area_EPA, 
     xlab = "PC1", 
     ylab = "PC2",
     col = cols(1000),
     xlim = c(-1.85, 2.88),
     ylim = c(-1.61, 7.17),
     cex.axis = 0.8,
     cex.lab = 0.8,
     axis.args = list(cex.axis = 0.8),
     font = 2, font.lab = 2
)
plot(pol, add = TRUE)

# writeRaster(env_space_area_EPA, "Processed_Data/env_space_EP.tif", format = "GTiff", overwrite = TRUE)

# Well-sampled cells ------------------------------------------------------

# WS_cent_2

# Import well-sampled cells with a slope threshold of 0.1 from the inventory completeness analysis
#well_sampled_2 <- read.csv("Data/well_sampled_2.csv")

# Select only well-sampled cells within the boundaries of Eastern Palearctic realm
well_sampled_EPA_2 <- well_sampled_2 %>% 
  dplyr::filter(Realm == "EP")

#_______________________________________________________________________________
# Now we calculate the environmental space for our Well-sampled cells
# Using estimators output from 'InventoryCompleteness' script:
well_sampled_EPA_2 <- well_sampled_EPA_2[, c('Longitude','Latitude')]
# Plot the well-sampled occurrences on the map
raster::plot(r_EPA, main = "", col = "gray", legend = FALSE)
points(well_sampled_EPA_2, col = rgb(1, 0, 0, .5), pch = 20, cex = 2)

# Extract climate values of WellSurvey cells
values_WS_EPA <- raster::extract(climate_E_Palearctic, 
                                 well_sampled_EPA_2,
                                 cellnumbers = TRUE)[, 1]

# write.table(values_WS_EPA, "Processed_Data/values_WS_EP.txt") # Save for environmental distance calculation

coords_WS_EPA <- climate_PCA_values_EPalearctic[values_WS_EPA, ]

# write.table(coords_WS_EPA, "Processed_Data/coords_WS_EP.txt") # Save for environmental distance calculation
cell_WS_EPA <- raster::extract(env_space_EPA, 
                               coords_WS_EPA,
                               cellnumbers = TRUE)[, 1]

n_WS_EPA <- table(cell_WS_EPA)
env_space_WS_EPA <- env_space_EPA
values(env_space_WS_EPA)[as.numeric(names(n_WS_EPA))] <- n_WS_EPA
ncell(cell_WS_EPA)

# The following map shows the frequency of climates in the well-sampled cells 
# in the study area (legend of this figure is continuous BUT values are integer numbers, 
# so manually adapt the legend to each case)

plot(env_space_WS_EPA, 
     xlab = "PC1", 
     ylab = "PC2",
     col = c("transparent", cols(1000)),
     xlim = c(-1.85, 2.88),
     ylim = c(-1.61, 7.17),
     font = 2, font.lab = 2,
     cex.axis = 0.8,
     cex.lab = 0.8,
     axis.args = list(cex.axis = 0.8))
plot(pol, add = TRUE)

# writeRaster(env_space_WS_EPA, "Processed_Data/env_space_WS_EP.tif", format = "GTiff", overwrite = TRUE)

# Schoener's D ------------------------------------------------------------

# Schoener's D: quantifies the overlap between the location of well-sampled 
# sites and cells with most frequent conditions 
# The Schoener's D index varies from zero (total lack of congruence) to one 
# (total congruence) 

# Transform the abundance of each cell into probabilities.
area_values_EPA <- area_values_EPA/sum(area_values_EPA, na.rm = TRUE) # Relative frequency of climate type for all the study area
WS_values_EPA <- values(env_space_WS_EPA)
WS_values_EPA <- WS_values_EPA/sum(WS_values_EPA, na.rm = TRUE) # Relative frequency of climate type for well-sampled cells

# Calculate the climate overlap using Schoener's D.
# D values close to 1 indicate that the location of well-sampled 
# sites coincide with climate conditions 
# that are frequently found in the study area
SchoenersD <- function(x, y) {
  sub_values <- abs(x - y)
  D <- 1 - (sum(sub_values, na.rm = TRUE) / 2)
  return(D)
}

D_EPA <- SchoenersD(area_values_EPA, WS_values_EPA)
print(paste("Climate overlap between well-sampled cells and the study area, given by the observed Schoener's D equals = ", 
            round(D_EPA,3), "%"))

# We create a null model to test if D values is different from a 
# random distribution of D values calculated from randomly sampling occurrence 
# records.
set.seed(0)
replications <- 1000 # Choose the number of replications
D_rnd <- numeric(replications)
for (i in 1:replications) {
  rnd <- sample(env_space_EPA_v, length(well_sampled_EPA_2), replace = TRUE)
  n_rnd <- table(rnd)
  env_space_rnd <- env_space_EPA
  values(env_space_rnd)[as.numeric(names(n_rnd))] <- n_rnd
  rnd_values <- values(env_space_rnd)
  rnd_values <- rnd_values/sum(rnd_values, na.rm = TRUE)
  D_rnd[i] <- SchoenersD(area_values_EPA, rnd_values)
}

# p-value from previous analysis
# If p < 0.05, it means that the location of well-sampled sites does not 
# coincide with areas with climate conditions frequently found in your study area
p_EPA <- (sum(D_EPA > D_rnd) + 1) / (length(D_rnd) + 1) # Unicaudal test
print(paste("p value equals = ", round(p_EPA, 3)))

# Plot distribution of model values
hist(D_rnd, 10, # Write number of bins
     xlim = c(0, 1), 
     main = (""),
     xlab = "D", 
     col = rgb(.5, .5, .5), 
     border = FALSE,
     font = 1, font.lab = 1)
abline(v = D_EPA, col = "red", lty = 2)# Red line show the observed Schoeners' D value

############ Kruskal-Wallis test ##################################
# The following map shows the distribution of well-sampled cells (red bars) 
# vs the study area cells (gray bars)
### For a) PC1
hist(myPCAEPalearctic$scores[, 1], 
     breaks = ncol(env_space_EPA),
     freq = F,
     col = "grey", 
     border = FALSE, 
     xlim= c(-3, 4), 
     ylim = c(0, 1.2), 
     main = "",
     xlab = "PC1", 
     font = 2, font.lab = 2, 
     cex.lab = 0.8, cex.axis = 0.8)
lines(density(na.omit(myPCAEPalearctic$scores[, 1])), col = "black", lwd = 2)

hist(coords_WS_EPA[, 1],
     breaks = 5,
     add = TRUE,
     col = adjustcolor("#0072B2", alpha.f = 0.5),
     freq = F,
     border = FALSE)
lines(density(na.omit(coords_WS_EPA[, 1])), col = "#0072B2", lwd = 2)

# X axis is a probability density
# Kruskal-Wallis verifies whether 1) the distribution of well-sampled sites 
# is an unbiased subset of the entire climate conditions of the Atlantic forest. 
# If this is so, p > 0.05 
x_1_EPA <- c(myPCAEPalearctic$scores[, 1], coords_WS_EPA[, 1])
g_1_EPA <- as.factor(c(rep("area", length(myPCAEPalearctic$scores[, 1])),
                       rep("WS", length(coords_WS_EPA[, 1]))))
kruskal.test(x_1_EPA ~ g_1_EPA)
# Kolmogorov smirnov test###
ks.test(myPCAEPalearctic$scores[, 1], coords_WS_EPA[, 1])

### For b) PC2
hist(myPCAEPalearctic$scores[, 2], 
     breaks = ncol(env_space_EPA),
     freq = F,
     col = "grey", 
     border = FALSE, 
     xlim = c(-2, 8), 
     ylim = c(0, 0.6), 
     main = "",
     xlab = "PC2", 
     font = 2, font.lab = 2, 
     cex.lab = 0.8, cex.axis = 0.8)
lines(density(na.omit(myPCAEPalearctic$scores[, 2])), col = "black", lwd = 2)

hist(coords_WS_EPA[, 2],
     breaks = 5, 
     add = TRUE,
     col = adjustcolor("#0072B2", alpha.f = 0.5),
     freq = F,
     border = FALSE)
lines(density(na.omit(coords_WS_EPA[, 2])), col = "#0072B2", lwd = 2)

x_2_EPA <- c(myPCAEPalearctic$scores[, 2], coords_WS_EPA[, 2])
g_2_EPA <- as.factor(c(rep("area", length(myPCAEPalearctic$scores[, 2])),
                       rep("WS", length(coords_WS_EPA[, 2]))))
kruskal.test(x_2_EPA ~ g_2_EPA)

# Kolmogorov smirnov test###
ks.test(myPCAEPalearctic$scores[, 2], coords_WS_EPA[, 2])

# Rarity ------------------------------------------------------------------

# We can check how many environmental space has been sampled 
# and how does these cells look like.
surface_EPA <- WS_values_EPA
surface_EPA[is.na(area_values_EPA) | area_values_EPA == 0] <- NA
sampled_EPA <- sum(surface_EPA > 0, na.rm = TRUE)/ sum(surface_EPA >= 0, na.rm = TRUE) * 100
percen_EPA <- round(sampled_EPA, 2)
prin_EPA <- paste0(percen_EPA, "% of our study area climate types covered by well-sampled cells")
print(prin_EPA)

# Is this env. space sampled corresponding to rare climates?
# First, we make values vary from 0 to 1 according to their rarity:
# Values close to 0, are very common, values close to 1 very rare
mini_EPA <- min(area_values_EPA, na.rm = TRUE) # less frequent value
# rarity index, also called Min-Max scalling
area_values01EPA <- abs(1 - (area_values_EPA - mini_EPA) / (max(area_values_EPA, na.rm = TRUE) - mini_EPA)) 
# see http://rasbt.github.io/mlxtend/user_guide/preprocessing/minmax_scaling/
hist(area_values01EPA, 
     breaks = ncol(env_space_EPA), 
     freq = FALSE, col = "white",
     main = "", 
     xlab = "Climate rarity", border = FALSE, 
     ylim = c(0, 10),
     font = 2, font.lab = 2,
     cex.lab = 0.8, cex.axis = 0.8)
lines(density(na.omit(area_values01EPA)), col = "black", lwd = 2)
lines(density(na.omit(area_values01EPA[surface_EPA > 0])), col = "#0072B2", lwd = 2)

# Kruskal-Wallis test to see if the distribution of rarities for the sampled 
# occurrences differs from the distribution observed in the entire area
x_EPA <- c(area_values01EPA, area_values01EPA[surface_EPA > 0])
g_EPA <- as.factor(c(rep("area", length(area_values01EPA)),
                     rep("ws", length(area_values01EPA[surface_EPA > 0]))))
kruskal.test(x_EPA ~ g_EPA)
# Kolmogorov smirnov test
area_values03_EPA <- (na.omit(area_values01EPA))
area_values04_EPA <- (na.omit(area_values01EPA[surface_EPA > 0]))
ks.test(area_values03_EPA, area_values04_EPA)

# This result indicates that the under-sampled area is composed 
# mainly by rare climates. You can inspect the environmental space figures 
# to check which areas were not sampled.

# Finally MAP the climatic rarity
rarity_env_EPA <- env_space_EPA
values(rarity_env_EPA) <- area_values01EPA
rarity_Percell_EPA <- raster::extract(rarity_env_EPA, climate_PCA_values_EPalearctic)
rarity_map_EPA <- r_EPA
values(rarity_map_EPA) <- rarity_Percell_EPA

cols <- colorRampPalette(c('#e0f3db','#a8ddb5','#43a2ca'))
x11()
plot(rarity_map_EPA, 
     col = cols(10), 
     font = 2, font.lab = 2, 
     ylim = c(23.2092, 81.2829), 
     xlim = c(-180, 180)
)

points(well_sampled_EPA_2, col = rgb(1, 0, 0, .5), pch = 19, cex = 1)

# Convert raster to data frame for ggplot2 if needed
rarity_map_df_EPA <- as.data.frame(as(rarity_map_EPA, "SpatialPixelsDataFrame"))
colnames(rarity_map_df_EPA) <- c("value", "x", "y")

# write.csv(rarity_map_df_EPA, "Processed_Data/Rarity_EP.csv") # Save to create rarity map
