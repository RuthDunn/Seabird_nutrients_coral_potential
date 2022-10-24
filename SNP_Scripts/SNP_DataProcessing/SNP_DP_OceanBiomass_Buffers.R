
rm(list = ls(all = TRUE))

# Packages:
library(raster)
library(rgdal)
library(dplyr)
library(ggplot2)
library(sf)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Read in ocean biomass data:

bio08 <- raster("Mapping/Ocean_biomass_mapping/Ocean_biomass_2008.tif")
bio09 <- raster("Mapping/Ocean_biomass_mapping/Ocean_biomass_2009.tif")
bio10 <- raster("Mapping/Ocean_biomass_mapping/Ocean_biomass_2010.tif")
bio11 <- raster("Mapping/Ocean_biomass_mapping/Ocean_biomass_2011.tif")
bio12 <- raster("Mapping/Ocean_biomass_mapping/Ocean_biomass_2012.tif")
bio13 <- raster("Mapping/Ocean_biomass_mapping/Ocean_biomass_2013.tif")
bio14 <- raster("Mapping/Ocean_biomass_mapping/Ocean_biomass_2014.tif")
bio15 <- raster("Mapping/Ocean_biomass_mapping/Ocean_biomass_2015.tif")
bio16 <- raster("Mapping/Ocean_biomass_mapping/Ocean_biomass_2016.tif")
bio17 <- raster("Mapping/Ocean_biomass_mapping/Ocean_biomass_2017.tif")
bio18 <- raster("Mapping/Ocean_biomass_mapping/Ocean_biomass_2018.tif")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Read in Chagos shapefile:

chagos <- readOGR("Mapping/Chagos_Shapefile.shp")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Transform both:

chagos <- spTransform(chagos, CRS("+proj=utm +zone=43 +datum=WGS84 +units=m +no_defs +type=crs"))
bio08 <- projectRaster(bio08, crs = crs(chagos))
bio09 <- projectRaster(bio09, crs = crs(chagos))
bio10 <- projectRaster(bio10, crs = crs(chagos))
bio11 <- projectRaster(bio11, crs = crs(chagos))
bio12 <- projectRaster(bio12, crs = crs(chagos))
bio13 <- projectRaster(bio13, crs = crs(chagos))
bio14 <- projectRaster(bio14, crs = crs(chagos))
bio15 <- projectRaster(bio15, crs = crs(chagos))
bio16 <- projectRaster(bio16, crs = crs(chagos))
bio17 <- projectRaster(bio17, crs = crs(chagos))
bio18 <- projectRaster(bio18, crs = crs(chagos))

# Check that this has worked okay by plotting it:

plot(bio08)
plot(chagos, add = T)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Create buffer areas: ####

buff300 <- buffer(chagos, width = 300000)
plot(buff300, add = T)

buff180 <- buffer(chagos, width = 180000)
plot(buff180, add = T)

buff600 <- buffer(chagos, width = 600000)
plot(buff600, add = T)

buff1200 <- buffer(chagos, width = 1200000)
plot(buff1200, add = T)

# Calculate areas
area300 <- area(buff300)
area180 <- area(buff180)
area600 <- area(buff600)
area1200 <- area(buff1200)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Extract raster values: ####

# 300 km buffer

vals300.bio08 <- extract(bio08, buff300, fun = sum)
vals300.bio09 <- extract(bio09, buff300, fun = sum)
vals300.bio10 <- extract(bio10, buff300, fun = sum)
vals300.bio11 <- extract(bio11, buff300, fun = sum)
vals300.bio12 <- extract(bio12, buff300, fun = sum)
vals300.bio13 <- extract(bio13, buff300, fun = sum)
vals300.bio14 <- extract(bio14, buff300, fun = sum)
vals300.bio15 <- extract(bio15, buff300, fun = sum)
vals300.bio16 <- extract(bio16, buff300, fun = sum)
vals300.bio17 <- extract(bio17, buff300, fun = sum)
vals300.bio18 <- extract(bio18, buff300, fun = sum)

mean300 <- mean(c(vals300.bio08, vals300.bio09, vals300.bio10, vals300.bio11, 
                  vals300.bio12, vals300.bio13, vals300.bio14, vals300.bio15, 
                  vals300.bio16, vals300.bio17, vals300.bio18))

sd300 <- sd(c(vals300.bio08, vals300.bio09, vals300.bio10, vals300.bio11, 
                  vals300.bio12, vals300.bio13, vals300.bio14, vals300.bio15, 
                  vals300.bio16, vals300.bio17, vals300.bio18))

rm(vals300.bio08, vals300.bio09, vals300.bio10, vals300.bio11, vals300.bio12, vals300.bio13, vals300.bio14,
   vals300.bio15, vals300.bio16, vals300.bio17, vals300.bio18)

# 180 km buffer

vals180.bio08 <- extract(bio08, buff180)
vals180.bio09 <- extract(bio09, buff180)
vals180.bio10 <- extract(bio10, buff180)
vals180.bio11 <- extract(bio11, buff180)
vals180.bio12 <- extract(bio12, buff180)
vals180.bio13 <- extract(bio13, buff180)
vals180.bio14 <- extract(bio14, buff180)
vals180.bio15 <- extract(bio15, buff180)
vals180.bio16 <- extract(bio16, buff180)
vals180.bio17 <- extract(bio17, buff180)
vals180.bio18 <- extract(bio18, buff180)

mean180 <- mean(c(sum(vals180.bio08[[1]]), sum(vals180.bio09[[1]]), sum(vals180.bio10[[1]]), sum(vals180.bio11[[1]]), 
                  sum(vals180.bio12[[1]]), sum(vals180.bio13[[1]]), sum(vals180.bio14[[1]]), sum(vals180.bio15[[1]]), 
                  sum(vals180.bio16[[1]]), sum(vals180.bio17[[1]]), sum(vals180.bio18[[1]])))

sd180 <- sd(c(sum(vals180.bio08[[1]]), sum(vals180.bio09[[1]]), sum(vals180.bio10[[1]]), sum(vals180.bio11[[1]]), 
              sum(vals180.bio12[[1]]), sum(vals180.bio13[[1]]), sum(vals180.bio14[[1]]), sum(vals180.bio15[[1]]), 
              sum(vals180.bio16[[1]]), sum(vals180.bio17[[1]]), sum(vals180.bio18[[1]])))

rm(vals180.bio08, vals180.bio09, vals180.bio10, vals180.bio11, vals180.bio12, vals180.bio13, vals180.bio14,
   vals180.bio15, vals180.bio16, vals180.bio17, vals180.bio18)

# 600 km buffer

vals600.bio08 <- extract(bio08, buff600)
vals600.bio09 <- extract(bio09, buff600)
vals600.bio10 <- extract(bio10, buff600)
vals600.bio11 <- extract(bio11, buff600)
vals600.bio12 <- extract(bio12, buff600)
vals600.bio13 <- extract(bio13, buff600)
vals600.bio14 <- extract(bio14, buff600)
vals600.bio15 <- extract(bio15, buff600)
vals600.bio16 <- extract(bio16, buff600)
vals600.bio17 <- extract(bio17, buff600)
vals600.bio18 <- extract(bio18, buff600)

mean600 <- mean(c(sum(vals600.bio08[[1]]), sum(vals600.bio09[[1]]), sum(vals600.bio10[[1]]), sum(vals600.bio11[[1]]), 
                  sum(vals600.bio12[[1]]), sum(vals600.bio13[[1]]), sum(vals600.bio14[[1]]), sum(vals600.bio15[[1]]), 
                  sum(vals600.bio16[[1]]), sum(vals600.bio17[[1]]), sum(vals600.bio18[[1]])))

sd600 <- sd(c(sum(vals600.bio08[[1]]), sum(vals600.bio09[[1]]), sum(vals600.bio10[[1]]), sum(vals600.bio11[[1]]), 
              sum(vals600.bio12[[1]]), sum(vals600.bio13[[1]]), sum(vals600.bio14[[1]]), sum(vals600.bio15[[1]]), 
              sum(vals600.bio16[[1]]), sum(vals600.bio17[[1]]), sum(vals600.bio18[[1]])))

rm(vals600.bio08, vals600.bio09, vals600.bio10, vals600.bio11, vals600.bio12, vals600.bio13, vals600.bio14,
   vals600.bio15, vals600.bio16, vals600.bio17, vals600.bio18)

# 1200 km buffer

vals1200.bio08 <- extract(bio08, buff1200)
vals1200.bio09 <- extract(bio09, buff1200)
vals1200.bio10 <- extract(bio10, buff1200)
vals1200.bio11 <- extract(bio11, buff1200)
vals1200.bio12 <- extract(bio12, buff1200)
vals1200.bio13 <- extract(bio13, buff1200)
vals1200.bio14 <- extract(bio14, buff1200)
vals1200.bio15 <- extract(bio15, buff1200)
vals1200.bio16 <- extract(bio16, buff1200)
vals1200.bio17 <- extract(bio17, buff1200)
vals1200.bio18 <- extract(bio18, buff1200)

mean1200 <- mean(c(sum(vals1200.bio08[[1]]), sum(vals1200.bio09[[1]]), sum(vals1200.bio10[[1]]), sum(vals1200.bio11[[1]]), 
                  sum(vals1200.bio12[[1]]), sum(vals1200.bio13[[1]]), sum(vals1200.bio14[[1]]), sum(vals1200.bio15[[1]]), 
                  sum(vals1200.bio16[[1]]), sum(vals1200.bio17[[1]]), sum(vals1200.bio18[[1]])))

sd1200 <- sd(c(sum(vals1200.bio08[[1]]), sum(vals1200.bio09[[1]]), sum(vals1200.bio10[[1]]), sum(vals1200.bio11[[1]]), 
              sum(vals1200.bio12[[1]]), sum(vals1200.bio13[[1]]), sum(vals1200.bio14[[1]]), sum(vals1200.bio15[[1]]), 
              sum(vals1200.bio16[[1]]), sum(vals1200.bio17[[1]]), sum(vals1200.bio18[[1]])))

rm(vals1200.bio08, vals1200.bio09, vals1200.bio10, vals1200.bio11, vals1200.bio12, vals1200.bio13, vals1200.bio14,
   vals1200.bio15, vals1200.bio16, vals1200.bio17, vals1200.bio18)

rm(bio08, bio09, bio10, bio11, bio12, bio13, bio14, bio15, bio16, bio17, bio18)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Tidy into dataframe: ####

df <- data.frame(mean = c(mean300, mean180, mean600, mean1200),
                 sd = c(sd300, sd180, sd600, sd1200),
                 area = c(area300, area180, area600, area1200))

rm(mean300, mean180, mean600, mean1200,
   sd300, sd180, sd600, sd1200,
   area300, area180, area600, area1200,
   buff300, buff180, buff600, buff1200)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# These values are in Jm-2
# So multiply them by area to find the total

df$meanj <- (df$area * df$mean)
df$sdj <- (df$area * df$sd)

df$meantonnes <- df$meanj/(4e+06)
df$sdtonnes <- df$sdj/(4e+06)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

write.csv(df, "SNP_Data/Processed/OceanBiomass_BufferAreas_R.csv")
