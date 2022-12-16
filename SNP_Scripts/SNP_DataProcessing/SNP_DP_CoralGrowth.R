rm(list = ls(all = TRUE))

library(plyr)
library("readxl")
library(tidyverse)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Load Casey's coral growth data ####

coral <- read_excel("SNP_Data/Raw data/Coral_Growth_Casey/Chagos_coral_growth_difference.xlsx", sheet = 1)

# Convert all character columns to factors

coral <- as.data.frame(unclass(coral),
                      stringsAsFactors = TRUE)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Load Nick's nitrogen input estimates

nitrogen.input <- read.csv("SNP_Data/Graham_etal_2018_Chagos_Seabirds_Nitrogen.csv")

# Convert all character columns to factors

nitrogen.input <- as.data.frame(unclass(nitrogen.input),
                                stringsAsFactors = TRUE)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Rename coral df islands

levels(coral$Island)[1] <- "Eagle_Island"
levels(coral$Island)[2] <- "Ile_Fouquet"
levels(coral$Island)[3] <- "Grande_Ile_Coquillage"
levels(coral$Island)[4] <- "Ile_Poule"
levels(coral$Island)[5] <- "Ile_Longue"
levels(coral$Island)[6] <- "Middle_Brother"
levels(coral$Island)[7] <- "PB_Ile_Anglaise"
levels(coral$Island)[8] <- "Sal_Ile_Anglaise"
levels(coral$Island)[9] <- "Ile_de_la_Passe"

# Combine dfs

coral <- join(coral, nitrogen.input[,c(2, 4, 11)], by = "Island") %>%
  select("Atoll", "Island", "diff_SA_year", "kg_N_ha_yr")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Write

write.csv(coral, "SNP_Data/Processed/CoralGrowth_Nitrogen.csv")
