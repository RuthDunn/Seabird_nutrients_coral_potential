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

# Load models ####

load("SNP_ModelOutputs/CoralGrowth_Nitrogen_c.Rdata")
load("SNP_ModelOutputs/Erosion_Nitrogen_c.Rdata")
load("SNP_ModelOutputs/Grazing_Nitrogen_c.Rdata")
load("SNP_ModelOutputs/Fish_Nitrogen_c.Rdata")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Load predication data ####

pred.data <- read.csv("SNP_Data/Processed/Seabird_NutrientInput_Predicted_OtherAtolls.csv")

# Calculate total nitrogen input (sum of the 3 species)
pred.data <- pred.data[,c("Atoll_Island", "Species", "NitrogenInput.perArea")]
# Long -> Wide
pred.data <- spread(pred.data, Species, NitrogenInput.perArea)
pred.data$NitrogenInput <- pred.data$Anous_tenuirostris + pred.data$Onychoprion_fuscatus + pred.data$Sula_sula

pred.data$logNitrogen <- log(pred.data$NitrogenInput)

# Center logNitrogen data around 1.86
pred.data$clogNitrogen <- pred.data$logNitrogen - 1.86

pred.data <- pred.data[,c("Atoll_Island", "logNitrogen")]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Make predictions ####

coral <- as.data.frame(predict(coralgrowth.model.run.center, newdata = pred.data, allow_new_levels = TRUE))
fish <- as.data.frame(predict(fishbiomass.model.run.center, newdata = pred.data, allow_new_levels = TRUE))
graz <- as.data.frame(predict(grazing.model.run.center, newdata = pred.data, allow_new_levels = TRUE))
eros <- as.data.frame(predict(erosion.model.run.center, newdata = pred.data, allow_new_levels = TRUE))

rm(coralgrowth.model.run.center, fishbiomass.model.run.center, grazing.model.run.center, erosion.model.run.center)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Combine into one dataframe

pred.data$coralgrowth <- coral$Estimate
pred.data$coralgrowth.low <- coral$Q2.5
pred.data$coralgrowth.high <- coral$Q97.5

pred.data$fishbiomass <- fish$Estimate
pred.data$fishbiomass.low <- fish$Q2.5
pred.data$fishbiomass.high <- fish$Q97.5

pred.data$grazing <- graz$Estimate
pred.data$grazing.low <- graz$Q2.5
pred.data$grazing.high <- graz$Q97.5

pred.data$erosion <- eros$Estimate
pred.data$erosion.low <- eros$Q2.5
pred.data$erosion.high <- eros$Q97.5

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

write.csv(pred.data, "SNP_Data/Processed/Pred_CoralFunction_OtherAtolls.csv")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Summary stats for ms

coral_max <- pred.data[which.max(pred.data$coralgrowth), ]
exp(coral_max[,c("coralgrowth", "coralgrowth.low", "coralgrowth.high")])

fish_max <- pred.data[which.max(pred.data$fishbiomass), ]
exp(fish_max[,c("fishbiomass", "fishbiomass.low", "fishbiomass.high")])

graz_max <- pred.data[which.max(pred.data$grazing), ]
exp(graz_max[,c("grazing", "grazing.low", "grazing.high")])

eros_max <- pred.data[which.max(pred.data$erosion), ]
exp(eros_max[,c("erosion", "erosion.low", "erosion.high")])
