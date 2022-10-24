rm(list = ls(all = TRUE))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Load data

veg.data <- read.csv("SNP_Data/Raw data/Tim_Wilkinson_RBGKew-Chagos_vegetation.csv")
# veg.categories <- read.csv("SNP_Data/Raw_data/Vegetation_Habitat_Wilkinson_Carr_Conversion.csv")

# Clean data

# Change NAs to 0s:
veg.data[is.na(veg.data)] <- 0

# Remove rows where the whole thing is 0s
veg.data <- veg.data[rowSums(veg.data[,4:21])>0,]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Edit data

# Create Savanna column:
veg.data$Savanna_ha <- veg.data$BareGround +
  veg.data$Grass + veg.data$Herb +
  veg.data$SandWithSparseCover

# Create Wetland column:
veg.data$Wetland_ha <- veg.data$BrackishWater +
  veg.data$Mangrove

# Create Native Forest column:
veg.data$NativeForest_ha <- veg.data$Broadleaf +
  veg.data$CordiaSp + veg.data$DeadCordia +
  veg.data$FernGrove + veg.data$MixedBroadleafCoconut +
  veg.data$PisoniaWoodland

# Create Non-Native Forest column:
veg.data$NonNativeForest_ha <- veg.data$CoconutPalm +
  veg.data$Papaya + veg.data$Unknown +
  veg.data$UnknownSp

# Create Mixed Shrub column:
veg.data$MixedShrub_ha <- veg.data$ScaevolaTaccada +
  veg.data$Thicket

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Select island, atoll, rat, and new habitat category cols only

veg.data <- veg.data[,c(1:3, 22:26)]

# How does Pete Carr use this data?
# % cover?
# Area (km2) and % cover
veg.data$Savanna_p <- veg.data$Savanna_ha/(veg.data$Savanna_ha + veg.data$Wetland_ha +
                                             veg.data$NativeForest_ha + veg.data$NonNativeForest_ha +
                                             veg.data$MixedShrub_ha) * 100

veg.data$Wetland_p <- veg.data$Wetland_ha/(veg.data$Savanna_ha + veg.data$Wetland_ha +
                                             veg.data$NativeForest_ha + veg.data$NonNativeForest_ha +
                                             veg.data$MixedShrub_ha) * 100

veg.data$NativeForest_p <- veg.data$NativeForest_ha/(veg.data$Savanna_ha + veg.data$Wetland_ha +
                                             veg.data$NativeForest_ha + veg.data$NonNativeForest_ha +
                                             veg.data$MixedShrub_ha) * 100

veg.data$NonNativeForest_p <- veg.data$NonNativeForest_ha/(veg.data$Savanna_ha + veg.data$Wetland_ha +
                                             veg.data$NativeForest_ha + veg.data$NonNativeForest_ha +
                                             veg.data$MixedShrub_ha) * 100

veg.data$MixedShrub_p <- veg.data$MixedShrub_ha/(veg.data$Savanna_ha + veg.data$Wetland_ha +
                                             veg.data$NativeForest_ha + veg.data$NonNativeForest_ha +
                                             veg.data$MixedShrub_ha) * 100

# Check that this is okay:
veg.data$Savanna_p + veg.data$Wetland_p + veg.data$NativeForest_p + veg.data$NonNativeForest_p + veg.data$MixedShrub_p
# Yep :)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

veg.data <- veg.data[complete.cases(veg.data),]

write.csv(veg.data, "SNP_Data/Habitat_cover.csv")
