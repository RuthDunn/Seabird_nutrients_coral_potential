rm(list = ls(all = TRUE))

library(tidyverse)

# Load raw data:

atoll.data <- read.csv("SNP_Data/Raw data/Islands_InvasiveSpecies_Atolls_Area.csv") %>%
  select("Atoll", "IslandArea")

# Condense into Atoll-specific values:

atoll.data <- as.data.frame(by(atoll.data$IslandArea,list(atoll.data$Atoll),sum))

# Write up:

write.csv(atoll.data, "SNP_Data/Processed/Tropical_Atolls_Invasives_Areas.csv")