rm(list = ls(all = TRUE))

library(tidyverse)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Load UVC data ####

fish <- read.csv("SNP_Data/Graham_etal_2018_Chagos_UVC_data.csv")

fish <- as.data.frame(unclass(fish),
                      stringsAsFactors = TRUE)

names(fish)[1] <- "Atoll"

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

nitrogen.input <- read.csv("SNP_Data/Graham_etal_2018_Chagos_Seabirds_Nitrogen.csv")

# Convert all character columns to factors

nitrogen.input <- as.data.frame(unclass(nitrogen.input),
                                stringsAsFactors = TRUE)

# Rename 'Eagle Island' as 'Eagle'

levels(nitrogen.input$Island)[1] <- "Eagle"

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Combine

fish <- join(fish, nitrogen.input[,c(2, 4, 11)], by = "Island") %>%
  select("Atoll", "Island", "Transect", "Biomass", "kg_N_ha_yr")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Write

write.csv(fish, "SNP_Data/Processed/FishBiomass_Nitrogen.csv")
