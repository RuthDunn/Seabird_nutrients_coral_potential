rm(list = ls(all = TRUE))

library(plyr)
library(ggplot2)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

seabirds <- read.csv("SNP_Data/Carr_etal_2020_Chagos_breeding_seabirds.csv")
names(seabirds)[1] <- "Island"

# Calc sum of all seabirds:
seabirds$Seabirds_Abundance_Total <- rowSums(seabirds[,c(4:21)])

# Remove rat islands:
# seabirds <- seabirds[which(seabirds$Rattus_rattus == "A"),]

# Atoll-Island column:
seabirds$Atoll_Island <- paste(seabirds$Atoll, seabirds$Island, sep = "_")

# Add in island areas:
areas <- read.csv("SNP_Data/Carr_2021_Chagos_islands_info.csv")
names(areas)[1] <- "Island"
areas$Atoll_Island <- paste(areas$Atoll, areas$Island, sep = "_")
seabirds <- join(seabirds, areas[,-c(1:3)], by = "Atoll_Island")
rm(areas)

# Add in habitat cover
habitat <- read.csv("SNP_Data/Wilkinson_2017_Chagos_habitat_cover.csv")
habitat <- habitat[,-1]
habitat$Atoll_Island <- paste(habitat$Atoll, "_", habitat[,1], sep = "")
seabirds <- join(seabirds, habitat[,-c(1:2)], by = "Atoll_Island")
rm(habitat)

seabirds <- seabirds[,c(23, 2,1,3, 24:36, 4:22)]

# Add in density metrics
# Total:
seabirds$Seabirds_Density_Total <- seabirds$Seabirds_Abundance_Total/seabirds$Size_Ha

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Write file
write.csv(seabirds, "SNP_Data/Processed/SeabirdDensityDrivers.csv")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
