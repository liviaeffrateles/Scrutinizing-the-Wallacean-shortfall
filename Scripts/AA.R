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
library(stringr)
library(openxlsx)

# Import the 19 CHELSA bioclimatic variables (1° resolution) cropped to the australasia realm

#Raster

path_env <- "D:/Scrutinizing-the-Wallacean-shortfall/Chelsa"
files_env <- list.files(path_env, pattern = "^climate_Australasia.*\\.tif$", full.names = TRUE)
allrasters <- lapply(files_env, raster)
clim_env <- raster::stack(files_env)

# Check the resulting aggregated raster to confirm the new resolution
res(clim_env)
r_AA <- raster(clim_env, 1) #extract 1 raster to check NEW resolution of cells
res_AA <- res(r_AA)  # Save value for later analysis
climate_Australasia <- clim_env # create climate stack of our new rasters
plot(climate_Australasia[[1]])

# PCA ---------------------------------------------------------------------

library(psych)
library(maps)
library(colorRamps)

# We reduce all the variables to fewer variables using a PCA.
# 1) Here we standardize and prepare the data.
valuesAustralasia <- as.data.frame(values(climate_Australasia)) # get env values from rasters

# 2) Standardization of the variables (Normalization - Mean = 0 and Std = 1)
std <- function(x){(x - mean(x, na.rm = T)) / sd(x, na.rm = T)} 
valuesAustralasiaStd <- apply(valuesAustralasia, 2, std) # apply to dataframe of climate values

rem <- apply(is.na(valuesAustralasiaStd), 1, any)
PCAAustralasia <- as.data.frame(valuesAustralasiaStd[!rem, ])

# 3) Run the PCA (2 axis)
matAustralasia <- matrix(runif(nrow(PCAAustralasia) * ncol(PCAAustralasia), 0.00001, 0.00009), 
                         ncol = ncol(PCAAustralasia)) # add a very small randomness to avoid singularity
PCAAustralasia2 <- PCAAustralasia + matAustralasia
myPCAAustralasia <- principal(PCAAustralasia2,
                              nfactors = 2,
                              rotate = "varimax",
                              scores = T)

# 4)
loadings_AA <- as.data.frame(unclass(myPCAAustralasia$loadings))

loadings_PC1_AA <- loadings_AA %>%
  mutate(abs_RC1 = abs(RC1)) %>%
  arrange(desc(abs_RC1))

loadings_PC2_AA <- loadings_AA %>%
  mutate(abs_RC2 = abs(RC2)) %>%
  arrange(desc(abs_RC2))

# 5) Plot PCA
biplot.psych(myPCAAustralasia, xlim.s = c(-2.2, 4), ylim.s = c(-4, 2)) # change limits to plot the whole pca values
rownames(loadings_AA) <- str_extract(rownames(loadings_AA), "bio\\d+")
loadings_AA <- tibble::rownames_to_column(loadings_AA, var = "variables")
scores_AA <- myPCAAustralasia$scores
scale_factor <- 2.5

biplot_AA <- ggplot() +
    geom_point(data = as.data.frame(scores_AA), aes(x = RC1, y = RC2), color = "black", size = 1) +
    geom_text(data = loadings_AA, aes(x = RC1  * scale_factor * 1.05, y = RC2  * scale_factor * 1.05, label = variables), 
              vjust = -0.5, hjust = 0.5, size = 2, color = "#F0A202", fontface = "bold") +
    geom_segment(data = loadings_AA, aes(x = 0, y = 0, xend = RC1 * scale_factor, yend = RC2 * scale_factor), 
                 arrow = arrow(type = "open", length = unit(0.1, "inches")), color = "#0072B2") +
    scale_x_continuous(breaks = seq(-2, 3, by = 1)) + 
    scale_y_continuous(breaks = seq(-4, 2, by = 1)) +
    xlab(paste("PC1 (56%)")) +
    ylab(paste("PC2 (24%)")) +
    theme_test(base_size = 12) +
    theme(aspect.ratio = 1)


biplot_AA

# PCA of our study area ---------------------------------------------------

# Instead of assuming the absolute values of all Chelsa variables, 
# the next map assumes values of the linear combinations between these variables (PCA scores).
valuesAustralasia_2 <- valuesAustralasiaStd[, 1:2] # create a vector with the same length as v2 but 2 columns
valuesAustralasia_2[!rem, ] <- myPCAAustralasia$scores 
climate_PCA_Australasia <- subset(climate_Australasia, 1:2)
values(climate_PCA_Australasia) <- valuesAustralasia_2 # insert pca values of v3 into raster
names(climate_PCA_Australasia) <- c("PC1", "PC2") # rename
# coords_PCA_AA <- as.data.frame(climate_PCA_Australasia, xy = T)
# write.table(coords_PCA_AA, "Processed_Data/coords_PCA_AA.txt") # Save for environmental distance calculation

# Plot study area by each axis from PCA values
world_sf <- ne_countries(scale = "medium", returnclass = "sf") %>%
  dplyr::filter(name != "Antarctica")

par(mfrow = c(1, 2))
plot(climate_PCA_Australasia, 1, font = 1, font.lab = 1, add = T)  # PCA axis 1
plot(world_sf, add = TRUE, border = "gray30", col = NA, lwd = 0.8, xlim = c(83, 180), ylim = c(-54.75, 4.8))
plot(climate_PCA_Australasia, 2, font = 1, font.lab = 1)  # PCA axis 2
plot(world_sf, add = TRUE, border = "gray30", col = NA, lwd = 0.8)

# Environmental space -----------------------------------------------------

# The next step is to create the environmental space using the two PCA scores.
# Create a Cartesian plan with PC1 and PC2 scores

climate_PCA_values_Australasia <- values(climate_PCA_Australasia) # get env values from PCA

# Transform the two vars we want into the env space, 
# by taking the min and max scores (PCA score values) for each PCA axis
# and create a raster object. 
xminAA <- min(climate_PCA_values_Australasia[, 1], na.rm = TRUE) 
xmaxAA <- max(climate_PCA_values_Australasia[, 1], na.rm = TRUE)
yminAA <- min(climate_PCA_values_Australasia[, 2], na.rm = TRUE)
ymaxAA <- max(climate_PCA_values_Australasia[, 2], na.rm = TRUE)

# This function creates the cartesian plan comprising the min and max PCA scores

env_space_AA <- raster(xmn = xminAA, xmx = xmaxAA, 
                       ymn = yminAA, ymx = ymaxAA, crs = 4326,
                       res = 0.5) # Resolution of the env cell can change

extent(env_space_AA) <- extent(climate_PCA_values_Australasia)
ncell(env_space_AA)
values(env_space_AA) <- 0
env_space_area_AA <- env_space_AA # duplicate this object for the next step

### Insert our PCA values into this env_space
# Extract PC1 and PC2, convert it in "classes of values" to plot in the map
env_space_AA_v <- raster::extract(env_space_AA, 
                                  y = na.omit(climate_PCA_values_Australasia), 
                                  cellnumbers = TRUE)[, 1]
table(env_space_AA_v)
n_env_space_AA_v <- table(env_space_AA_v) # this function counts the frequency of "climates" (PC scores)


values(env_space_area_AA)[as.numeric(names(n_env_space_AA_v))] <- n_env_space_AA_v
area_values_AA <- values(env_space_area_AA)
area_values_AA[area_values_AA == 0] <- NA
values(env_space_area_AA) <- area_values_AA
pol <- rasterToPolygons(env_space_area_AA, dissolve = FALSE)

# The following map shows the frequency of climates in the study area.
# The axes shows all the PC1 and PC2 scores. 
# Each cell "roughly" represents a climate type.
# The colors indicate the frequency of climate types.

cols <- colorRampPalette(c("#56B4E9", "#33A02C", "#FDE725FF", "#CC79A7", "#D55E00"))

plot(env_space_area_AA, 
       xlab = "PC1", 
       ylab = "PC2",
       col = cols(1000),
       xlim = c(-1.14, 3.5),
       ylim = c(-4.52, 1.82),
       font = 2, font.lab = 2,
       cex.axis = 0.8, cex.lab = 0.8,
       axis.args = list(cex.axis = 0.8))
  
plot(pol, add = TRUE)

# writeRaster(env_space_area_AA, "Processed_Data/env_space_AA.tif", format = "GTiff", overwrite = TRUE)

# Well-sampled cells ------------------------------------------------------

# Import well-sampled cells with a slope threshold of 0.05 from the inventory completeness analysis
# well_sampled_1 <- read.csv("Processed_Data/well_sampled_1.csv")

# Select only well-sampled cells within the boundaries of Australasia realm
well_sampled_AA_1 <- well_sampled_1 %>% 
  dplyr::filter(Realm == "AA")

#_______________________________________________________________________________
# Now we calculate the environmental space for our Well-sampled cells
# Using estimators output from 'InventoryCompleteness' script:
well_sampled_AA_1 <- well_sampled_AA_1[, c('Longitude','Latitude')]
# Plot the well-sampled occurrences on the map
raster::plot(r_AA, main = "", col = "gray", legend = FALSE)
points(well_sampled_AA_1, col = rgb(1, 0, 0, .5), pch = 20, cex = 2)

# Extract climate values of Well-sampled cells
values_WS_AA <- raster::extract(climate_Australasia, 
                                well_sampled_AA_1,
                                cellnumbers = TRUE)[, 1]
# write.table(values_WS_AA, "Processed_Data/values_WS_AA.txt") # Save for environmental distance calculation

coords_WS_AA <- climate_PCA_values_Australasia[values_WS_AA, ]
# write.table(coords_WS_AA, "Processed_Data/coords_WS_AA.txt") # Save for environmental distance calculation

cell_WS_AA <- raster::extract(env_space_AA, 
                              coords_WS_AA,
                              cellnumbers = TRUE)[, 1]

n_WS_AA <- table(cell_WS_AA)
env_space_WS_AA <- env_space_AA
values(env_space_WS_AA)[as.numeric(names(n_WS_AA))] <- n_WS_AA
ncell(cell_WS_AA)

# The following map shows the frequency of climates in the well-sampled cells 
# in the study area (legend of this figure is continuous BUT values are integer numbers, 
# so manually adapt the legend to each case)

plot(env_space_WS_AA, 
     xlab = "PC1", 
     ylab = "PC2",
     xlim = c(-1.14, 3.5),
     ylim = c(-4.52, 1.82),
     col = c("transparent", cols(1000)),
     font = 2, font.lab = 2,
     cex.axis = 0.8, cex.lab = 0.8,
     axis.args = list(cex.axis = 0.8))
plot(pol, add = TRUE)

# writeRaster(env_space_WS_AA, "Processed_Data/env_space_WS_AA.tif", format = "GTiff", overwrite = TRUE)

# Schoener's D ------------------------------------------------------------

# Schoener's D: quantifies the overlap between the location of well-sampled 
# sites and cells with most frequent conditions 
# The Schoener's D index varies from zero (total lack of congruence) to one 
# (total congruence) 

# Transform the abundance of each cell into probabilities.
area_values_AA <- area_values_AA/sum(area_values_AA, na.rm = TRUE) # Relative frequency of climate type for all the study area
WS_values_AA <- values(env_space_WS_AA)
WS_values_AA <- WS_values_AA/sum(WS_values_AA, na.rm = TRUE) # Relative frequency of climate type for well-sampled cells

# Calculate the climate overlap using Schoener's D.
# D values close to 1 indicate that the location of well-sampled 
# sites coincide with climate conditions 
# that are frequently found in the study area
SchoenersD <- function(x, y) {
  sub_values <- abs(x - y)
  D <- 1 - (sum(sub_values, na.rm = TRUE) / 2)
  return(D)
}

D_AA <- SchoenersD(area_values_AA, WS_values_AA)
print(paste("Climate overlap between well-sampled cells and the study area, given by the observed Schoener's D equals = ", 
            round(D_AA, 3), "%"))

# We create a null model to test if D values is different from a 
# random distribution of D values calculated from randomly sampling occurrence 
# records.
set.seed(0)
replications <- 1000 # Choose the number of replications
D_rnd <- numeric(replications)
for (i in 1:replications) {
  rnd <- sample(env_space_AA_v, length(well_sampled_AA_1), replace = TRUE)
  n_rnd <- table(rnd)
  env_space_rnd <- env_space_AA
  values(env_space_rnd)[as.numeric(names(n_rnd))] <- n_rnd
  rnd_values <- values(env_space_rnd)
  rnd_values <- rnd_values/sum(rnd_values, na.rm = TRUE)
  D_rnd[i] <- SchoenersD(area_values_AA, rnd_values)
}

# p-value from previous analysis
# If p < 0.05, it means that the location of well-sampled sites does not 
# coincide with areas with climate conditions frequently found in your study area
p_AA <- (sum(D_AA > D_rnd) + 1) / (length(D_rnd) + 1) # Unicaudal test
print(paste("p value equals = ", round(p_AA, 3)))

# Plot distribution of model values
hist(D_rnd, 10, # Write number of bins
     xlim = c(0, 1), 
     main = (""),
     xlab = "D", 
     col = rgb(.5, .5, .5), 
     border = FALSE,
     font = 1, font.lab = 1)
abline(v = D_AA, col = "red", lty = 2) # Red line show the observed Schoeners' D value

# Kruskal-Wallis test -----------------------------------------------------

# The following map shows the distribution of well-sampled cells (red bars) 
# vs the study area cells (gray bars)
### For a) PC1

hist(myPCAAustralasia$scores[, 1], 
     breaks = ncol(env_space_AA),
     freq = F,
     col = "grey", 
     border = FALSE, 
     xlim = c(-1.78, 4.10), 
     ylim = c(0, 1.45), 
     main = "",
     xlab = "PC1", 
     font = 2, font.lab = 2, 
     cex.lab = 0.8, cex.axis = 0.8)
lines(density(na.omit(myPCAAustralasia$scores[, 1])), col = "black", lwd = 2)

hist(coords_WS_AA[, 1],
     add = TRUE,
     col = adjustcolor("#0072B2", alpha.f = 0.5),
     freq = F,
     border = FALSE)
lines(density(na.omit(coords_WS_AA[, 1])), col = "#0072B2", lwd = 2)

# X axis is a probability density
# Kruskal-Wallis verifies whether 1) the distribution of well-sampled sites 
# is an unbiased subset of the entire climate conditions of the Neotropical Region. 
# If this is so, p > 0.05 
x_1_AA <- c(myPCAAustralasia$scores[, 1], coords_WS_AA[, 1])
g_1_AA <- as.factor(c(rep("area", length(myPCAAustralasia$scores[, 1])),
                      rep("WS", length(coords_WS_AA[, 1]))))
kruskal.test(x_1_AA ~ g_1_AA)
# Kolmogorov smirnov test###
ks.test(myPCAAustralasia$scores[, 1], coords_WS_AA[, 1])

### For b) PC2
hist(myPCAAustralasia$scores[, 2], 
     breaks = ncol(env_space_AA),
     freq = F,
     col = "grey", 
     border = FALSE, 
     xlim = c(-5.7, 3), 
     ylim = c(0, 0.5), 
     main = "",
     xlab = "PC2", 
     font = 2, font.lab = 2, 
     cex.lab = 0.8, cex.axis = 0.8)
lines(density(na.omit(myPCAAustralasia$scores[, 2])), col = "black", lwd = 2)

hist(coords_WS_AA[, 2],
     breaks = 5, 
     add = TRUE,
     col = adjustcolor("#0072B2", alpha.f = 0.5),
     freq = F,
     border = FALSE)
lines(density(na.omit(coords_WS_AA[, 2])), col = "#0072B2", lwd = 2)

x_2_AA <- c(myPCAAustralasia$scores[, 2], coords_WS_AA[, 2])
g_2_AA <- as.factor(c(rep("area", length(myPCAAustralasia$scores[, 2])),
                      rep("WS", length(coords_WS_AA[, 2]))))

kruskal.test(x_2_AA ~ g_2_AA)

# Kolmogorov smirnov test###
ks.test(myPCAAustralasia$scores[, 2], coords_WS_AA[, 2])

# Rarity ------------------------------------------------------------------

# We can check how many environmental space has been sampled 
# and how does these cells look like.
surface_AA <- WS_values_AA
surface_AA[is.na(area_values_AA) | area_values_AA == 0] <- NA
sampled_AA <- sum(surface_AA > 0, na.rm = TRUE)/ sum(surface_AA >= 0, na.rm = TRUE) * 100
percen_AA <- round(sampled_AA, 2)
prin_AA <- paste0(percen_AA, "% of our study area climate types covered by well-sampled cells")
print(prin_AA)

# Is this env. space sampled corresponding to rare climates?
# First, we make values vary from 0 to 1 according to their rarity:
# Values close to 0, are very common, values close to 1 very rare
mini_AA <- min(area_values_AA, na.rm = TRUE) # less frequent value
# rarity index, also called Min-Max scalling
area_values01AA <- abs(1 - (area_values_AA - mini_AA) / (max(area_values_AA, na.rm = TRUE) - mini_AA)) 
# see http://rasbt.github.io/mlxtend/user_guide/preprocessing/minmax_scaling/

hist(area_values01AA, 
     breaks = ncol(env_space_AA), 
     freq = FALSE, col = "white",
     main = "", 
     xlab = "Climate rarity", border = FALSE, 
     ylim = c(0, 5),
     font = 2, font.lab = 2,
     cex.lab = 0.8, cex.axis = 0.8)
lines(density(na.omit(area_values01AA)), col = "black", lwd = 2)
lines(density(na.omit(area_values01AA[surface_AA > 0])), col = "#0072B2", lwd = 2)

# Kruskal-Wallis test to see if the distribution of rarities for the sampled 
# occurrences differs from the distribution observed in the entire area
x_AA <- c(area_values01AA, area_values01AA[surface_AA > 0])
g_AA <- as.factor(c(rep("area", length(area_values01AA)),
                    rep("ws", length(area_values01AA[surface_AA > 0]))))
kruskal.test(x_AA ~ g_AA)
# Kolmogorov smirnov test
area_values03_AA <- (na.omit(area_values01AA))
area_values04_AA <- (na.omit(area_values01AA[surface_AA > 0]))
ks.test(area_values03_AA, area_values04_AA)

# This result indicates that the under-sampled area is composed 
# mainly by rare climates. You can inspect the environmental space figures 
# to check which areas were not sampled.

# Finally MAP the climatic rarity
rarity_env_AA <- env_space_AA
values(rarity_env_AA) <- area_values01AA
rarity_Percell_AA <- raster::extract(rarity_env_AA, climate_PCA_values_Australasia)
rarity_map_AA <- r_AA
values(rarity_map_AA) <- rarity_Percell_AA

cols <- colorRampPalette(c('#e0f3db','#a8ddb5','#43a2ca'))
x11()
plot(rarity_map_AA, 
     col = cols(10), 
     font = 2, font.lab = 2, 
     ylim = c(-54.75187, 4.776925), 
     xlim = c(100, 179.202)
)
points(well_sampled_AA_1, col = "red", pch = 19, cex = 0.7)

# Convert raster to data frame 
rarity_map_df_AA <- as.data.frame(as(rarity_map_AA, "SpatialPixelsDataFrame"))
colnames(rarity_map_df_AA) <- c("value", "x", "y")

# write.csv(rarity_map_df_AA, "Processed_Data/Rarity_AA.csv") # Save to create rarity maps







