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

# Import the 19 CHELSA bioclimatic variables (1° resolution) cropped to the Neotropical realm

#Raster
path_env <- "D:/Scrutinizing-the-Wallacean-shortfall/Chelsa"
files_env <- list.files(path_env, pattern = "^climate_Neotropical.*\\.tif$", full.names = TRUE)
allrasters <- lapply(files_env, raster)
clim_env <- raster::stack(files_env)

# Check the resulting aggregated raster to confirm the new resolution
res(clim_env)
r_NT <- raster(clim_env, 1) #extract 1 raster to check NEW resolution of cells
res_NT <- res(r_NT)  # Save value for later analysis
climate_Neotropical <- clim_env # create climate stack of our new rasters
plot(climate_Neotropical[[1]])

# PCA ---------------------------------------------------------------------

library(psych)
library(maps)
library(colorRamps)

# We reduce all the variables to fewer variables using a PCA.
# 1) Here we standardize and prepare the data.
valuesNeotropical <- as.data.frame(values(climate_Neotropical)) # get env values from rasters

# 2) Standardization of the variables (Normalization - Mean = 0 and Std = 1)
std <- function(x){(x - mean(x, na.rm = T)) / sd(x, na.rm = T)} 
valuesNeotropicalStd <- apply(valuesNeotropical, 2, std) # apply to dataframe of climate values

rem <- apply(is.na(valuesNeotropicalStd), 1, any)
PCANeotropical <- as.data.frame(valuesNeotropicalStd[!rem, ])

# 3) Run the PCA (2 axis)
matNeotropical <- matrix(runif(nrow(PCANeotropical) * ncol(PCANeotropical), 0.00001, 0.00009), 
                         ncol = ncol(PCANeotropical)) # add a very small randomness to avoid singularity
PCANeotropical2 <- PCANeotropical + matNeotropical
myPCANeotropical <- principal(PCANeotropical2,
                              nfactors = 2,
                              rotate = "varimax",
                              scores = T)

# 4)

loadings_NT <- as.data.frame(unclass(myPCANeotropical$loadings))

loadings_PC1_NT <- loadings_NT%>%
  mutate(abs_RC1 = abs(RC1)) %>%
  arrange(desc(abs_RC1))

loadings_PC2_NT <- loadings_NT %>%
  mutate(abs_RC2 = abs(RC2)) %>%
  arrange(desc(abs_RC2))

# 5) Plot PCA
biplot.psych(myPCANeotropical, xlim.s = c(-4, 2), ylim.s = c(-2, 7)) # change limits to plot the whole pca values
rownames(loadings_NT) <- str_extract(rownames(loadings_NT), "bio\\d+")
loadings_NT <- tibble::rownames_to_column(loadings_NT, var = "variables")
scores_NT <- myPCANeotropical$scores
scale_factor <- 2.5

biplot_NT <- ggplot() +
    geom_point(data = as.data.frame(scores_NT), aes(x = RC1, y = RC2), color = "black", size = 1) +
    geom_segment(data = loadings_NT, aes(x = 0, y = 0, xend = RC1 * scale_factor, yend = RC2 * scale_factor), 
                 arrow = arrow(type = "open", length = unit(0.1, "inches")), color = "#0072B2") +
    geom_text(data = loadings_NT, aes(x = RC1  * scale_factor * 1.05, y = RC2  * scale_factor * 1.05, label = variables), 
              vjust = -0.5, hjust = 0.5, size = 2, color = "#F0A202", fontface = "bold") +
    scale_x_continuous(breaks = seq(-4, 2, by = 2)) + 
    scale_y_continuous(breaks = seq(-2, 7, by = 2)) +
    xlab(paste("PC1 (40%)")) +
    ylab(paste("PC2 (31%)")) +
    theme_test(base_size = 12) +
    theme(aspect.ratio = 1)

biplot_NT

# PCA of our study area ---------------------------------------------------

# Instead of assuming the absolute values of all Chelsa variables, 
# the next map assumes values of the linear combinations between these variables (PCA scores).
valuesNeotropical_2 <- valuesNeotropicalStd[, 1:2] # create a vector with the same length as v2 but 2 columns
valuesNeotropical_2[!rem, ] <- myPCANeotropical$scores 
climate_PCA_Neotropical <- subset(climate_Neotropical, 1:2)
values(climate_PCA_Neotropical) <- valuesNeotropical_2 # insert pca values of v3 into raster
names(climate_PCA_Neotropical) <- c("PC1", "PC2") # rename
# coords_PCA_NT <- as.data.frame(climate_PCA_Neotropical, xy = T) 
# write.table(coords_PCA_NT, "Processed_Data/coords_PCA_NT.txt") # Save for environmental distance calculation

# Plot study area by each axis from PCA values
world_sf <- ne_countries(scale = "medium", returnclass = "sf") %>%
  dplyr::filter(name != "Antarctica")

par(mfrow = c(1, 2))
plot(climate_PCA_Neotropical, 1, font = 1, font.lab = 1, add = T)  # PCA axis 1
plot(world_sf, add = TRUE, border = "gray30", col = NA, lwd = 0.8, xlim = c(-138, -8), ylim = c(-58, 33))
plot(climate_PCA_Neotropical, 2, font = 1, font.lab = 1)  # PCA axis 2
plot(world_sf, add = TRUE, border = "gray30", col = NA, lwd = 0.8)

# Environmental space -----------------------------------------------------

# The next step is to create the environmental space using the two PCA scores.
# Create a Cartesian plan with PC1 and PC2 scores

climate_PCA_values_Neotropical <- values(climate_PCA_Neotropical) # get env values from PCA

# Transform the two vars we want into the env space, 
# by taking the min and max scores (PCA score values) for each PCA axis
# and create a raster object. 
xminNT <- min(climate_PCA_values_Neotropical[, 1], na.rm = TRUE) 
xmaxNT <- max(climate_PCA_values_Neotropical[, 1], na.rm = TRUE)
yminNT <- min(climate_PCA_values_Neotropical[, 2], na.rm = TRUE)
ymaxNT <- max(climate_PCA_values_Neotropical[, 2], na.rm = TRUE)

# This function creates the cartesian plan comprising the min and max PCA scores

env_space_NT <- raster(xmn = xminNT, xmx = xmaxNT, 
                       ymn = yminNT, ymx = ymaxNT, crs = 4326,
                       res = 0.5) # Resolution of the env cell can change

extent(env_space_NT) <- extent(climate_PCA_values_Neotropical)

ncell(env_space_NT)
values(env_space_NT) <- 0
env_space_area_NT <- env_space_NT # duplicate this object for the next step

# Insert our PCA values into this env_space
# Extract PC1 and PC2, convert it in "classes of values" to plot in the map
env_space_NT_v <- raster::extract(env_space_NT, 
                                  y = na.omit(climate_PCA_values_Neotropical), 
                                  cellnumbers = TRUE)[, 1]
table(env_space_NT_v)
n_env_space_NT_v <- table(env_space_NT_v) # this function counts the frequency of "climates" (PC scores)

values(env_space_area_NT)[as.numeric(names(n_env_space_NT_v))] <- n_env_space_NT_v
area_values_NT <- values(env_space_area_NT)
area_values_NT[area_values_NT == 0] <- NA
values(env_space_area_NT) <- area_values_NT
pol <- rasterToPolygons(env_space_area_NT, dissolve = FALSE)

# The following map shows the frequency of climates in the study area.
# The axes shows all the PC1 and PC2 scores. 
# Each cell "roughly" represents a climate type.
# The colors indicate the frequency of climate types.

cols <- colorRampPalette(c("#56B4E9", "#33A02C", "#FDE725FF", "#CC79A7", "#D55E00"))

plot(env_space_area_NT, 
     xlab = "PC1", 
     ylab = "PC2",
     col = cols(1000),
     xlim = c(-4, 1.3),
     ylim = c(-1.7, 6.9),
     font = 2, font.lab = 2,
     cex.axis = 0.8,
     cex.lab = 0.8,
     axis.args = list(cex.axis = 0.8)
)
plot(pol, add = TRUE)

# writeRaster(env_space_area_NT, "Processed_Data/env_space_NT.tif", format = "GTiff", overwrite = TRUE)

# Well-sampled cells ------------------------------------------------------

# Import well-sampled cells with a slope threshold of 0.1 from the inventory completeness analysis
# well_sampled_2 <- read.csv("Data/well_sampled_2.csv")

# Select only well-sampled cells within the boundaries of Neotropical realm
well_sampled_NT_2 <- well_sampled_2 %>% 
  dplyr::filter(Realm == "NT")

#_______________________________________________________________________________
# Now we calculate the environmental space for our well-sampled cells cells
# Using estimators output from 'InventoryCompleteness' script:
well_sampled_NT_2 <- well_sampled_NT_2[, c('Longitude','Latitude')]
# Plot the well-sampled occurrences on the map
raster::plot(r_NT, main = "", col = "gray", legend = FALSE)
points(well_sampled_NT_2, col = rgb(1, 0, 0, .5), pch = 20, cex = 2)

# Extract climate values of Well-sampled cells
values_WS_NT <- raster::extract(climate_Neotropical, 
                                well_sampled_NT_2,
                                cellnumbers = TRUE)[, 1]

# write.table(values_WS_NT, "Processed_Data/values_WS_NT.txt") # Save for environmental distance calculation

coords_WS_NT <- climate_PCA_values_Neotropical[values_WS_NT, ]

# write.table(coords_WS_NT, "Processed_Data/coords_WS_NT.txt") # Save for environmental distance calculation
cell_WS_NT <- raster::extract(env_space_NT, 
                              coords_WS_NT,
                              cellnumbers = TRUE)[, 1]

n_WS_NT <- table(cell_WS_NT)
env_space_WS_NT <- env_space_NT
values(env_space_WS_NT)[as.numeric(names(n_WS_NT))] <- n_WS_NT
ncell(cell_WS_NT)

# The following map shows the frequency of climates in the well-sampled cells 
# in the study area (legend of this figure is continuous BUT values are integer numbers, 
# so manually adapt the legend to each case)
plot(env_space_WS_NT,
     xlab = "PC1",
     ylab = "PC2",
     col = c("transparent", cols(1000)),
     xlim = c(-4, 1.3),
     ylim = c(-1.7, 6.9),
     font = 2, font.lab = 2,
     cex.axis = 0.8, cex.lab = 0.8,
     axis.args = list(cex.axis = 0.8))
plot(pol, add = TRUE)

# writeRaster(env_space_WS_NT, "Processed_Data/env_space_WS_NT.tif", format = "GTiff", overwrite = TRUE)

# Schoener's D ------------------------------------------------------------

# Schoener's D: quantifies the overlap between the location of well-sampled 
# sites and cells with most frequent conditions 
# The Schoener's D index varies from zero (total lack of congruence) to one 
# (total congruence) 

# Transform the abundance of each cell into probabilities.
area_values_NT <- area_values_NT/sum(area_values_NT, na.rm = TRUE) # Relative frequency of climate type for all the study area
WS_values_NT <- values(env_space_WS_NT)
WS_values_NT <- WS_values_NT/sum(WS_values_NT, na.rm = TRUE) # Relative frequency of climate type for well-sampled cells

# Calculate the climate overlap using Schoener's D.
# D values close to 1 indicate that the location of well-sampled 
# sites coincide with climate conditions 
# that are frequently found in the study area
SchoenersD <- function(x, y) {
  sub_values <- abs(x - y)
  D <- 1 - (sum(sub_values, na.rm = TRUE) / 2)
  return(D)
}

D_NT <- SchoenersD(area_values_NT, WS_values_NT)
print(paste("Climate overlap between well-sampled cells and the study area, given by the observed Schoener's D equals = ", 
            round(D_NT, 3), "%"))

# We create a null model to test if D values is different from a 
# random distribution of D values calculated from randomly sampling occurrence 
# records.
set.seed(0)
replications <- 1000 # Choose the number of replications
D_rnd <- numeric(replications)
for (i in 1:replications) {
  rnd <- sample(env_space_NT_v, length(well_sampled_NT_2), replace = TRUE)
  n_rnd <- table(rnd)
  env_space_rnd <- env_space_NT
  values(env_space_rnd)[as.numeric(names(n_rnd))] <- n_rnd
  rnd_values <- values(env_space_rnd)
  rnd_values <- rnd_values/sum(rnd_values, na.rm = TRUE)
  D_rnd[i] <- SchoenersD(area_values_NT, rnd_values)
}

# p-value from previous analysis
# If p < 0.05, it means that the location of well-sampled sites does not 
# coincide with areas with climate conditions frequently found in your study area
p_NT <- (sum(D_NT > D_rnd) + 1) / (length(D_rnd) + 1) # Unicaudal test
print(paste("p value equals = ", round(p_NT, 3)))

# Plot distribution of model values
hist(D_rnd, 10, # Write number of bins
     xlim = c(0, 1), 
     main = (""),
     xlab = "D", 
     col = rgb(.5, .5, .5), 
     border = FALSE,
     font = 1, font.lab = 1)
abline(v = D_NT, col = "red", lty = 2) # Red line show the observed Schoeners' D value

# Kruskal-Wallis test -----------------------------------------------------

# The following map shows the distribution of well-sampled cells (red bars) 
# vs the study area cells (gray bars)
### For a) PC1
hist(myPCANeotropical$scores[, 1], 
     breaks = ncol(env_space_NT),
     freq = F,
     col = "grey", 
     border = FALSE, 
     xlim= c(-4, 2), 
     ylim = c(0, 0.8), 
     main = "",
     xlab = "PC1", 
     font = 2, font.lab = 2, 
     cex.lab = 0.8, cex.axis = 0.8)
lines(density(na.omit(myPCANeotropical$scores[, 1])), col = "black", lwd = 2)

hist(coords_WS_NT[, 1],
     add = TRUE,
     col = adjustcolor("#0072B2", alpha.f = 0.5),
     freq = F,
     border = FALSE)
lines(density(na.omit(coords_WS_NT[, 1])), col = "#0072B2", lwd = 2)

# X axis is a probability density
# Kruskal-Wallis verifies whether 1) the distribution of well-sampled sites 
# is an unbiased subset of the entire climate conditions of the Neotropical Region. 
# If this is so, p > 0.05 
x_1_NT <- c(myPCANeotropical$scores[, 1], coords_WS_NT[, 1])
g_1_NT <- as.factor(c(rep("area", length(myPCANeotropical$scores[, 1])),
                      rep("WS", length(coords_WS_NT[, 1]))))
kruskal.test(x_1_NT ~ g_1_NT)
# Kolmogorov smirnov test###
ks.test(myPCANeotropical$scores[, 1], coords_WS_NT[, 1])

### For b) PC2
hist(myPCANeotropical$scores[, 2], 
     breaks = ncol(env_space_NT),
     freq = F,
     col = "grey", 
     border = FALSE, 
     xlim = c(-2, 8), 
     ylim = c(0, 0.6), 
     main = "",
     xlab = "PC2", 
     font = 2, font.lab = 2, 
     cex.lab = 0.8, cex.axis = 0.8)
lines(density(na.omit(myPCANeotropical$scores[, 2])), col = "black", lwd = 2)

hist(coords_WS_NT[, 2],
     breaks = 5, 
     add = TRUE,
     col = adjustcolor("#0072B2", alpha.f = 0.5),
     freq = F,
     border = FALSE)
lines(density(na.omit(coords_WS_NT[, 2])), col = "#0072B2", lwd = 2)

x_2_NT <- c(myPCANeotropical$scores[, 2], coords_WS_NT[, 2])
g_2_NT <- as.factor(c(rep("area", length(myPCANeotropical$scores[, 2])),
                      rep("WS", length(coords_WS_NT[, 2]))))

kruskal.test(x_2_NT ~ g_2_NT)

# Kolmogorov smirnov test###
ks.test(myPCANeotropical$scores[, 2], coords_WS_NT[, 2])

# Rarity ------------------------------------------------------------------

# We can check how many environmental space has been sampled 
# and how does these cells look like.
surface_NT <- WS_values_NT
surface_NT[is.na(area_values_NT) | area_values_NT == 0] <- NA
sampled_NT <- sum(surface_NT > 0, na.rm = TRUE)/ sum(surface_NT >= 0, na.rm = TRUE) * 100
percen_NT <- round(sampled_NT, 2)
prin_NT <- paste0(percen_NT, "% of our study area climate types covered by well-sampled cells")
print(prin_NT)

# Is this env. space sampled corresponding to rare climates?
# First, we make values vary from 0 to 1 according to their rarity:
# Values close to 0, are very common, values close to 1 very rare
mini_NT <- min(area_values_NT, na.rm = TRUE) # less frequent value
# rarity index, also called Min-Max scalling
area_values01NT <- abs(1 - (area_values_NT - mini_NT) / (max(area_values_NT, na.rm = TRUE) - mini_NT)) 
# see http://rasbt.github.io/mlxtend/user_guide/preprocessing/minmax_scaling/
hist(area_values01NT, 
     breaks = ncol(env_space_NT), 
     freq = FALSE, col = "white",
     main = "", 
     xlab = "Climate rarity", border = FALSE, 
     ylim = c(0, 12),
     font = 2, font.lab = 2,
     cex.lab = 0.8, cex.axis = 0.8)
lines(density(na.omit(area_values01NT)), col = "black", lwd = 2)
lines(density(na.omit(area_values01NT[surface_NT > 0])), col = "#0072B2", lwd = 2)

# Kruskal-Wallis test to see if the distribution of rarities for the sampled 
# occurrences differs from the distribution observed in the entire area
x_NT <- c(area_values01NT, area_values01NT[surface_NT > 0])
g_NT <- as.factor(c(rep("area", length(area_values01NT)),
                    rep("ws", length(area_values01NT[surface_NT > 0]))))
kruskal.test(x_NT ~ g_NT)
# Kolmogorov smirnov test
area_values03_NT <- (na.omit(area_values01NT))
area_values04_NT <- (na.omit(area_values01NT[surface_NT > 0]))
ks.test(area_values03_NT, area_values04_NT)

# This result indicates that the under-sampled area is composed 
# mainly by rare climates. You can inspect the environmental space figures 
# to check which areas were not sampled.

# Finally MAP the climatic rarity

rarity_env_NT <- env_space_NT
values(rarity_env_NT) <- area_values01NT
rarity_Percell_NT <- raster::extract(rarity_env_NT, climate_PCA_values_Neotropical)
rarity_map_NT <- r_NT
values(rarity_map_NT) <- rarity_Percell_NT


cols <- colorRampPalette(c('#e0f3db','#a8ddb5','#43a2ca'))
x11()
plot(rarity_map_NT, 
     col = cols(10), 
     font = 2, font.lab = 2, 
     ylim = c(-55.97751, 29.36736), 
     xlim = c(-114.8031, -29.27657)
)
points(well_sampled_NT_2, col = "red", pch = 19, cex = 0.7)

# Convert raster to data frame 
rarity_map_df_NT <- as.data.frame(as(rarity_map_NT, "SpatialPixelsDataFrame"))
colnames(rarity_map_df_NT) <- c("value", "x", "y")

# write.csv(rarity_map_df_NT, "Processed_Data/Rarity_NT.csv") # Save to create rarity map

