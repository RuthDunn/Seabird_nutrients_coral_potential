rm(list = ls(all = TRUE))

library(ggplot2)
library(RColorBrewer)
library(patchwork)

library(sf)
library(ggspatial)
library(rnaturalearth)

library(scales)
library(plyr)
library(ggpubr)
library(tidyr)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Chagos centroid ####

# Load:
chagos <- st_read("Mapping/Chagos_Centroid.shp")

# Convert coordinates to "WGS 84 / PDC Mercator":
chagos <- st_transform(chagos, crs = 3832)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Atoll centroids ####

# Load:
centroids <- st_read("Mapping/Selected_Atoll_Invasive_Islands/Selected_Atoll_Centroids.shp")

# Convert coordinates to "WGS 84 / PDC Mercator":
centroids <- st_transform(centroids, crs = 3832)

centroids.df <- as.data.frame(centroids)

# st_crs(centroids)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Join up with prediction data ####

# A) Seabirds
seabirds <- read.csv("SNP_Data/Processed/Seabird_NutrientInput_Predicted_OtherAtolls.csv")
seabirds$Atoll <- seabirds$Atoll_Island

# Calculate total nitrogen input (sum of the 3 species)
seabirds <- seabirds[,c("Atoll_Island", "Species", "Estimate", "tot.Nitrogen")]
# Long -> Wide
birds <- spread(seabirds[,c("Atoll_Island", "Species", "Estimate")], Species, Estimate)
nitrogen <- spread(seabirds[,c("Atoll_Island", "Species", "tot.Nitrogen")], Species, tot.Nitrogen)

# B) Nitrogen
coralfunc <- read.csv("SNP_Data/Processed/Pred_CoralFunction_OtherAtolls.csv")

# Combine:
coralfunc$seabirds <- birds$Anous_tenuirostris + birds$Onychoprion_fuscatus + birds$Sula_sula
coralfunc$TotalNitrogen <- nitrogen$Anous_tenuirostris + nitrogen$Onychoprion_fuscatus + nitrogen$Sula_sula
coralfunc$Atoll <- coralfunc$Atoll_Island

rm(seabirds, birds, nitrogen)

# Convert to proportional increase:
coralfunc$coralgrowth <- rescale(coralfunc$coralgrowth)
coralfunc$fishbiomass <- rescale(coralfunc$fishbiomass)
coralfunc$grazing <- rescale(coralfunc$grazing)
coralfunc$erosion <- rescale(coralfunc$erosion)
# Calculate mean of 4 metrics instead:
coralfunc$mean.coral.metric <- rowMeans(coralfunc[,c("coralgrowth", "fishbiomass", "grazing", "erosion")])

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Transform back into spatial object:

atoll.predictions <- join(coralfunc[,-c(1:2)],
                          centroids.df[,c("Atoll", "geometry")], by = "Atoll")

atoll.predictions <- st_as_sf(atoll.predictions)

rm(coralfunc, centroids.df)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# World map ####

# Load:
worldmap <- ne_countries(returnclass = "sf")

# Reproject:
worldmap <- st_transform(worldmap, crs = 3857)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Seabirds & Nitrogen Plot ####

(summary(atoll.predictions$TotalNitrogen))

p1 <- ggplot() +
  geom_hline(yintercept = 0, col = "grey70", linetype = "dashed") +
  geom_sf(data = worldmap, fill = "#009E73", col = "#009E73") +
  geom_sf(data = chagos, col = "#FFB000", size = 4, shape = 8) +
  
  # Prediction data
  geom_sf(data = atoll.predictions,
          aes(colour = TotalNitrogen, size = seabirds), alpha = 0.7) +
  
  # Make pretty
  scale_colour_gradient2(low = "#648FFF", mid = "#DC267F", high = "#FFB000", midpoint = 15000,
                         name = expression("Seabird-derived nitrogen input (kg year"^-1*")"),
                         # name = expression("Seabird-derived nitrogen input (kg year"^-1*"ha"^-1*")"),
                         labels = scales::comma) +
  scale_size(labels = scales::comma, name = "Seabird abundance (number of pairs)") +
  # scale_alpha(range = c(0.4, 1)) +
  coord_sf(crs = st_crs(3832), xlim = c(-17000000 , 17000000 ), ylim = c(-5000000 , 5000000)) +
  guides(colour = guide_colourbar(barwidth = 15, barheight = 1)) +
  theme_light() %+replace% theme(legend.position="bottom", legend.box="vertical", legend.margin=margin(),
        axis.text.x = element_text(),
        axis.text.y = element_text(),
        plot.margin = margin(0, 0, 0, 0, unit = "pt")) +
  ylab(element_blank())

p1

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Coral Functions Plot ####

p2 <- ggplot() +
  geom_hline(yintercept = 0, col = "grey70", linetype = "dashed") +
  geom_sf(data = worldmap, fill = "#009E73", col = "#009E73") +
  geom_sf(data = chagos, col = "#FFB000", size = 4, shape = 8) +
  
  # Prediction data
  geom_sf(data = atoll.predictions,
          aes(colour = mean.coral.metric), size = 3, alpha = 0.7) +

  # Make pretty
  # Make pretty
  scale_colour_gradient2(low = "#648FFF", mid = "#DC267F", high = "#FFB000", midpoint = 0.5,
                         name = expression("Mean relative benefit to coral reef metrics")) +
  # scale_alpha(range = c(0.4, 1)) +
  coord_sf(crs = st_crs(3832), xlim = c(-17000000 , 17000000 ), ylim = c(-5000000 , 5000000)) +
  guides(colour = guide_colourbar(barwidth = 15, barheight = 1)) +
  theme_light() %+replace% theme(legend.position="bottom", legend.box="vertical", legend.margin=margin(),
                                 axis.text.x = element_text(),
                                 axis.text.y = element_text(),
                                 plot.margin = margin(0, 0, 0, 0, unit = "pt")) +
  ylab(element_blank())

p2

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ggarrange(p1, p2, nrow = 2, heights = c(3,2.5))

ggsave("Plots/Map_Atolls_Predictions_Everything.png", width = 8, height = 7)

ggsave("Plots/Map_Atolls_Predictions_Everything_poster.png", width = 40, height = 30, units = "cm")
