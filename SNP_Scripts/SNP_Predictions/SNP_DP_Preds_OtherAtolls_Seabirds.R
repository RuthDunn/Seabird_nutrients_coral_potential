rm(list = ls(all = TRUE))

# ~~~~~~~~~~~~~~~~~~~~~~~

# Packages:
library(brms)
library(dplyr)
library(tidybayes)
library(ggplot2)
library(patchwork)
library(ggmcmc)
library(jtools)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Load other atoll rat data ####

islands <- read.csv("SNP_Data/Processed/Tropical_Atolls_Invasives_Areas.csv")
head(islands)

# ~~~~~~~~~~~~~~~~~~~~~~~

# Load Chagos model ####

load("SNP_ModelOutputs/Seabirds_Habitat_brms_priors.Rdata")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Modify the data so that we can make predictions ####

# Requirements:
# Atoll_Island; Rattus_rattus; logArea; Species; id; logVeg

islands$Atoll_Island <- islands$Atoll
islands$Rattus_rattus <- "A"
# Convert area from km2 to ha and log
islands$logArea <- log(islands$Area*100)

islands.long <- do.call("rbind", replicate(3,
                                           islands[,c("Atoll_Island", "Rattus_rattus", "logArea")],
                                           simplify = FALSE))

islands.long$Species <- rep(c("Sula_sula", "Onychoprion_fuscatus", "Anous_tenuirostris"), each = nrow(islands))
islands.long$id <- as.numeric(rownames(islands.long))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Make predictions based on model parameters ####

# 1) If we eradicate rats but keep high proportions of plantation forest
# (the lower quartile of non plantation forest cover in Chagos: 26%)

islands.long$NativeVeg <- 25

predictions.highnn <- exp(fitted(seabirds.model, newdata = islands.long, allow_new_levels = TRUE, dpar = "mu"))

predictions.highnn <- cbind(islands.long,predictions.highnn)

# 2) If we eradicate rats and have low proportions of plantation forest
# (the UPPER quartile of NON-plantation forest cover in Chagos: 73%)

islands.long$NativeVeg <- 75

predictions.lownn <- exp(fitted(seabirds.model, newdata = islands.long, allow_new_levels = TRUE, dpar = "mu"))

predictions.lownn <- cbind(islands.long,predictions.lownn)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ggplot(predictions.lownn, aes(x = Atoll_Island, y = (Estimate), ymin = (Q2.5), ymax = (Q97.5))) +
  geom_point() +
  geom_errorbar() +
  facet_grid(Species~.) +
  scale_y_continuous(labels = scales::comma)

ggplot(predictions.lownn, aes(x = Estimate)) + geom_histogram() + facet_grid(Species~.)
ggplot(predictions.lownn, aes(x = Q97.5)) + geom_histogram() + facet_grid(Species~.)
ggplot(predictions.lownn, aes(x = Q2.5)) + geom_histogram() + facet_grid(Species~.)

# Densities
predictions.lownn$Estimate/exp(predictions.lownn$logArea)

# Save up:

write.csv(predictions.highnn, "SNP_ModelOutputs/OtherAtolls_Seabird_Predictions_HighNNForest.hln.csv")
write.csv(predictions.lownn, "SNP_ModelOutputs/OtherAtolls_Seabird_Predictions_LowNNForest.hln.csv")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Summary stats for ms

birds_max <- predictions.lownn[which.max(predictions.lownn$Estimate), ]
birds_max

area_small <- predictions.lownn[which.min(predictions.lownn$logArea), ]
area_small
