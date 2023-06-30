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

pred.data <- read.csv("SNP_Data/Processed/Seabird_NutrientInput_Chagos_Predicted_priors.csv")

# Calculate total nitrogen input (sum of the 3 species)
pred.data <- pred.data[,c("Atoll_Island", "Species", "NitrogenInput.Current", "NitrogenInput.high.nn", "NitrogenInput.low.nn")]

# Find total nitrogen input per Atoll_Island
# (Under the different scenarios)
pred.data <- aggregate(pred.data[,c(3:5)], by = list(pred.data$Atoll_Island), FUN = "sum")

# Remove DG
pred.data <- pred.data[-1,]

# Rename Atoll_Island col
colnames(pred.data)[1] <- "Atoll_Island"

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Get data ready to make predictions:

# Current predictions:

pred.data$logNitrogen <- log(pred.data$NitrogenInput.Current)

# Center logNitrogen data around 1.86
pred.data$clogNitrogen <- pred.data$logNitrogen - 1.86

pred.data.current <- pred.data[,c("Atoll_Island", "clogNitrogen")]

coral.current <- as.data.frame(predict(coralgrowth.model.run.center, newdata = pred.data.current, allow_new_levels = TRUE))
fish.current <- as.data.frame(predict(fishbiomass.model.run.center, newdata = pred.data.current, allow_new_levels = TRUE))
graz.current <- as.data.frame(predict(grazing.model.run.center, newdata = pred.data.current, allow_new_levels = TRUE))
eros.current <- as.data.frame(predict(erosion.model.run.center, newdata = pred.data.current, allow_new_levels = TRUE))

# Bad veg predictions:

pred.data$logNitrogen <- log(pred.data$NitrogenInput.high.nn)

# Center logNitrogen data around 1.86
pred.data$clogNitrogen <- pred.data$logNitrogen - 1.86

pred.data.bveg <- pred.data[,c("Atoll_Island", "clogNitrogen")]

coral.bveg <- as.data.frame(predict(coralgrowth.model.run.center, newdata = pred.data.bveg, allow_new_levels = TRUE))
fish.bveg <- as.data.frame(predict(fishbiomass.model.run.center, newdata = pred.data.bveg, allow_new_levels = TRUE))
graz.bveg <- as.data.frame(predict(grazing.model.run.center, newdata = pred.data.bveg, allow_new_levels = TRUE))
eros.bveg <- as.data.frame(predict(erosion.model.run.center, newdata = pred.data.bveg, allow_new_levels = TRUE))

# Good veg predictions:

pred.data$logNitrogen <- log(pred.data$NitrogenInput.low.nn)

# Center logNitrogen data around 1.86
pred.data$clogNitrogen <- pred.data$logNitrogen - 1.86

pred.data.gveg <- pred.data[,c("Atoll_Island", "clogNitrogen")]

coral.gveg <- as.data.frame(predict(coralgrowth.model.run.center, newdata = pred.data.gveg, allow_new_levels = TRUE))
fish.gveg <- as.data.frame(predict(fishbiomass.model.run.center, newdata = pred.data.gveg, allow_new_levels = TRUE))
graz.gveg <- as.data.frame(predict(grazing.model.run.center, newdata = pred.data.gveg, allow_new_levels = TRUE))
eros.gveg <- as.data.frame(predict(erosion.model.run.center, newdata = pred.data.gveg, allow_new_levels = TRUE))

rm(coralgrowth.model.run.center, erosion.model.run.center, grazing.model.run.center, fishbiomass.model.run.center)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Combine into one dataframe

pred.data <- pred.data[,-5]

pred.data$coralgrowth.current <- coral.current$Estimate
pred.data$coralgrowth.current.low <- coral.current$Q2.5
pred.data$coralgrowth.current.high <- coral.current$Q97.5

pred.data$coralgrowth.bveg <- coral.bveg$Estimate
pred.data$coralgrowth.bveg.low <- coral.bveg$Q2.5
pred.data$coralgrowth.bveg.high <- coral.bveg$Q97.5

pred.data$coralgrowth.gveg <- coral.gveg$Estimate
pred.data$coralgrowth.gveg.low <- coral.gveg$Q2.5
pred.data$coralgrowth.gveg.high <- coral.gveg$Q97.5

#

pred.data$fishbiomass.current <- fish.current$Estimate
pred.data$fishbiomass.current.low <- fish.current$Q2.5
pred.data$fishbiomass.current.high <- fish.current$Q97.5

pred.data$fishbiomass.bveg <- fish.bveg$Estimate
pred.data$fishbiomass.bveg.low <- fish.bveg$Q2.5
pred.data$fishbiomass.bveg.high <- fish.bveg$Q97.5

pred.data$fishbiomass.gveg <- fish.gveg$Estimate
pred.data$fishbiomass.gveg.low <- fish.gveg$Q2.5
pred.data$fishbiomass.gveg.high <- fish.gveg$Q97.5

#

pred.data$grazing.current <- graz.current$Estimate
pred.data$grazing.current.low <- graz.current$Q2.5
pred.data$grazing.current.high <- graz.current$Q97.5

pred.data$grazing.bveg <- graz.bveg$Estimate
pred.data$grazing.bveg.low <- graz.bveg$Q2.5
pred.data$grazing.bveg.high <- graz.bveg$Q97.5

pred.data$grazing.gveg <- graz.gveg$Estimate
pred.data$grazing.gveg.low <- graz.gveg$Q2.5
pred.data$grazing.gveg.high <- graz.gveg$Q97.5

#

pred.data$erosion.current <- eros.current$Estimate
pred.data$erosion.current.low <- eros.current$Q2.5
pred.data$erosion.current.high <- eros.current$Q97.5

pred.data$erosion.bveg <- eros.bveg$Estimate
pred.data$erosion.bveg.low <- eros.bveg$Q2.5
pred.data$erosion.bveg.high <- eros.bveg$Q97.5

pred.data$erosion.gveg <- eros.gveg$Estimate
pred.data$erosion.gveg.low <- eros.gveg$Q2.5
pred.data$erosion.gveg.high <- eros.gveg$Q97.5

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Write:

write.csv(pred.data, "SNP_Data/Processed/Prediction_data/Pred_CoralFunction_Chagos_c.csv")
