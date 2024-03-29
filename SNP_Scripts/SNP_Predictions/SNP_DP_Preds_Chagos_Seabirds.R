rm(list = ls(all = TRUE))

# Packages:
library(brms)
library(dplyr)
library(tidybayes)
library(ggplot2)
library(patchwork)
library(ggmcmc)
library(jtools)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Load Chagos rat data and model ####

islands <- read.csv("SNP_Data/Processed/ChagosRatIslands.csv")
islands <- islands[,-1]

# Remove DG main island:

islands <- islands[-1,]

# Load model of the influences of habitat on seabird abundance:

load("SNP_ModelOutputs/Seabirds_Habitat_brms_priors.Rdata")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Modify the data so that we can make predictions ####

islands$Rattus_rattus <- "A"
islands$logArea <- log(islands$Size_Ha)

islands.long <- do.call("rbind", replicate(3, islands, simplify = FALSE))

islands.long$Species <- rep(c("Sula_sula", "Onychoprion_fuscatus", "Anous_tenuirostris"), each = nrow(islands))
islands.long$id <- as.numeric(rownames(islands.long))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Make predictions based on model parameters ####

# 1) If we eradicate rats but keep high proportions of plantation forest
# (the lower quartile of non plantation forest cover in Chagos: 26%)
islands.long$NativeVeg <- 25

predictions1 <- exp(fitted(seabirds.model, newdata = islands.long[,c("Atoll_Island", "Species",
                                                                         "Rattus_rattus", "logArea", "NativeVeg")],
                           allow_new_levels = TRUE, dpar = "mu"))

predictions1 <- cbind(islands.long,predictions1)

# 2) If we eradicate rats but keep high proportions of plantation forest
# (the lower quartile of non plantation forest cover in Chagos: 73%)
islands.long$NativeVeg <- 50

predictions2 <- exp(fitted(seabirds.model, newdata = islands.long[,c("Atoll_Island", "Species",
                                                                     "Rattus_rattus", "logArea", "NativeVeg")],
                           allow_new_levels = TRUE, dpar = "mu"))

predictions2 <- cbind(islands.long,predictions2)

# 3) If we eradicate rats but keep high proportions of plantation forest
# (the lower quartile of non plantation forest cover in Chagos: 73%)
islands.long$NativeVeg <- 75

predictions3 <- exp(fitted(seabirds.model, newdata = islands.long[,c("Atoll_Island", "Species",
                                                                         "Rattus_rattus", "logArea", "NativeVeg")],
                           allow_new_levels = TRUE, dpar = "mu"))

predictions3 <- cbind(islands.long,predictions3)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

write.csv(predictions1, "SNP_ModelOutputs/Chagos_Seabird_Predictions_HighNNForest_priors.csv")
write.csv(predictions2, "SNP_ModelOutputs/Chagos_Seabird_Predictions_MidNNForest_priors.csv")
write.csv(predictions3, "SNP_ModelOutputs/Chagos_Seabird_Predictions_LowNNForest_priors.csv")
