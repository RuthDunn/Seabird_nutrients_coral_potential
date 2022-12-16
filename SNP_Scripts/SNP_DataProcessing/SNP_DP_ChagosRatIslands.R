rm(list = ls(all = TRUE))

library(plyr)
library(ggplot2)
library(tidyverse)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Load island area and habitat cover data

# Island areas:
areas <- read.csv("SNP_Data/Carr_2021_Chagos_islands_info.csv")
names(areas)[1] <- "Island"
areas$Atoll_Island <- paste(areas$Atoll, areas$Island, sep = "_")

# Habitat cover:
habitat <- read.csv("SNP_Data/Wilkinson_2017_Chagos_habitat_cover.csv")
habitat <- habitat[,-1]
habitat$Atoll_Island <- paste(habitat$Atoll, "_", habitat[,1], sep = "")

# Join:
island.data <- join(areas, habitat[,c(13,11)], by = "Atoll_Island")
rm(areas, habitat)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Keep ratty islands only

levels(as.factor(island.data$Rats))

island.data <- island.data[which(island.data$Rats=='P'),] %>%
  select("Island", "Atoll", "Rats", "Size_Ha", "Atoll_Island", "NonNativeForest_p")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Write file

write.csv(island.data, "SNP_Data/Processed/ChagosRatIslands.csv")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
