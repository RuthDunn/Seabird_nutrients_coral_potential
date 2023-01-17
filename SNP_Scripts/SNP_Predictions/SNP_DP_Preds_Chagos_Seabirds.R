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

load("SNP_ModelOutputs/Seabirds_Habitat_brms.hln.Rdata")

seabirds.model.run <- seabirds.model.run11
rm(seabirds.model.run11)

# pp_check(seabirds.model.run)

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
islands.long$logVeg <- log(25)

predictions1 <- exp(fitted(seabirds.model.run, newdata = islands.long[,c("Atoll_Island", "Species", "Size_Ha",
                                                                         "Rattus_rattus", "logArea", "logVeg")],
                           allow_new_levels = TRUE, dpar = "mu"))

predictions1 <- cbind(islands.long,predictions1)

# 2) If we eradicate rats but keep high proportions of plantation forest
# (the lower quartile of non plantation forest cover in Chagos: 73%)
islands.long$logVeg <- log(75)

predictions2 <- exp(fitted(seabirds.model.run, newdata = islands.long[,c("Atoll_Island", "Species", "Size_Ha",
                                                                         "Rattus_rattus", "logArea", "logVeg")],
                           allow_new_levels = TRUE, dpar = "mu"))

predictions2 <- cbind(islands.long,predictions2)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ggplot(predictions1, aes(x = Atoll_Island, y = (Estimate), ymin = (Q2.5), ymax = (Q97.5))) +
  geom_point() +
  geom_errorbar() +
  facet_grid(Species~.) +
  scale_y_continuous(labels = scales::comma)

ggplot(predictions1, aes(x = Estimate)) + geom_histogram() + facet_grid(Species~.)
ggplot(predictions1, aes(x = Q97.5)) + geom_histogram() + facet_grid(Species~.)
ggplot(predictions1, aes(x = Q2.5)) + geom_histogram() + facet_grid(Species~.)

ggplot(predictions2, aes(x = Atoll_Island, y = (Estimate), ymin = (Q2.5), ymax = (Q97.5))) +
  geom_point() +
  geom_errorbar() +
  facet_grid(Species~.) +
  scale_y_continuous(labels = scales::comma)

ggplot(predictions2, aes(x = Estimate)) + geom_histogram() + facet_grid(Species~.)
ggplot(predictions2, aes(x = Q97.5)) + geom_histogram() + facet_grid(Species~.)
ggplot(predictions2, aes(x = Q2.5)) + geom_histogram() + facet_grid(Species~.)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

write.csv(predictions1, "SNP_ModelOutputs/Chagos_Seabird_Predictions_HighNNForest.hln.csv")
write.csv(predictions2, "SNP_ModelOutputs/Chagos_Seabird_Predictions_LowNNForest.hln.csv")
