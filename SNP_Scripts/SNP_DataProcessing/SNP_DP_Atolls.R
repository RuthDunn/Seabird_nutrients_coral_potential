rm(list = ls(all = TRUE))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Load data

atoll.data <- read.csv("SNP_Data/Raw data/Islands_InvasiveSpecies_Atolls_Area.csv")
atoll.data <- atoll.data[,c("Atoll", "IslandArea", "IslandCoas")]
atoll.data["Atoll"] <- lapply(atoll.data["Atoll"] , factor)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Condense into Atoll-specific values

Area <- (by(atoll.data$IslandArea,list(atoll.data$Atoll),sum))
Coastline <- (by(atoll.data$IslandCoas,list(atoll.data$Atoll),sum))

new <- cbind(Area, Coastline)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

write.csv(new, "SNP_Data/Processed/Tropical_Atolls_Invasives_Areas.csv")
