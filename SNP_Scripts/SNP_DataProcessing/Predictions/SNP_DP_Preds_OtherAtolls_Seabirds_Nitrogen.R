rm(list = ls(all = TRUE))

library(plyr)
library(ggplot2)
library(ape)
library(MCMCglmm)
library(data.table)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Data loads ####

# Load predicted seabird numbers:
preds <- read.csv("SNP_ModelOutputs/OtherAtolls_Seabird_Predictions_LowNNForest.hln.csv")
head(preds)

# Long to wide:
# preds <- gather(preds[,-c(1,7)], Species, Estimate, Anous_tenuirostris:Sula_sula)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Nitrogen input per ha (NI) calculation:

# (Graham et al 2018)
# NI[i,j] = N[g] * Dr[i] * Bd[i,j] * Res[i,j] / IsArea[j]

# N[g] = nitrogen content of guano = 18.1% 
# Dr[i] = daily defecation rate (g) per species i ~ scaled allometrically from boobies
# Bd = number of birds of species i on each island j
#       - Sooty tern = 1 adult per pair present each 24 hour period
#       - Rf booby = 1 adult per pair present each 24 hour period + 1 adult present for 12 hours
#       - Lesser noddy = 2 adults for each 24 hour period
# Res = number of days of the year that the species is resident
# IsArea = island area

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Bring in:
#       - Island area
#       - Seabird breeding season length
#       - Seabird mass

# 1) Island area:
areas <- read.csv("SNP_Data/Processed/Tropical_Atolls_Invasives_Areas.csv")
areas$Atoll_Island <- areas$Atoll
# convert from km2 to ha
areas$Area <- areas$Area_km2 * 100
preds <- join(preds, areas[,c("Atoll_Island", "Area")], by = "Atoll_Island")
tot.area.ha <- colSums(areas[7])
rm(areas)


# 2) Seabird breeding season length
# 3) Seabird mass
bird.info <- read.csv("SNP_Data/Schreiber_Burger_2002_Seabird_info.csv")
names(bird.info)[1] <- "Species"
bird.info$Season.length <- bird.info$IncPhase + bird.info$BroodPhase + bird.info$CrechePhase
seabirds <- join(preds, bird.info[,c("Species", "Season.length", "Mass")], by = "Species")
rm(bird.info, preds)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Calculate defecation rate
# "We found birds defecated approximately every 36 min, each dropping weighing 0.66 g"
# So daily defecation rate = 
# (1440/36) * 0.66 = 26.4 g/day
# "For calculating defecation quantities for other species, we assume guano produced scaled
# allometrically with body size of species." (Young et al., 2010)
# DR[i] = DR[booby] * Mass[i]^0.63/Mass[booby]^0.63 (Staunton Smith & Johnson 1995)

seabirds$Defecation.rate <- (26.4 * (seabirds$Mass*1000)^0.63)/1150^0.63
# grams per day

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# NI[i,j] = N[g] * Dr[i] * Bd[i,j] * Res[i,j] / IsArea[j]
# Ng = nitrogen content of guano
# Dr = defecation rate (in g) per species (i) of bird per day
# Bd = number of that species of bird (i) on the island (j)
# Res = the number of days that the bird is resident
# IsArea = island area:

# Adjust Bd (number of days on the island) to account for time away whilst foraging
# Sooty tern - one bird foraging - assume one adult present
# Red-footed booby - one adult always present, one away over night - assume 1.5 birds present
# Lesser noddy - short foraging trips, don't adjust numbers - assume 2 birds present
seabirds$Estimate.adjusted <- ifelse(seabirds$Species == "Sula_sula", seabirds$Estimate * 1.5,
                                     ifelse(seabirds$Species == "Anous_tenuirostris", seabirds$Estimate * 2,
                            seabirds$Estimate))
seabirds$Estimate.adjusted.low <- ifelse(seabirds$Species == "Sula_sula", seabirds$Q2.5 * 1.5,
                                     ifelse(seabirds$Species == "Anous_tenuirostris", seabirds$Q2.5 * 2,
                                            seabirds$Q2.5))
seabirds$Estimate.adjusted.high <- ifelse(seabirds$Species == "Sula_sula", seabirds$Q97.5 * 1.5,
                                     ifelse(seabirds$Species == "Anous_tenuirostris", seabirds$Q97.5 * 2,
                                            seabirds$Q97.5))

# (Adjust defecation rate in g to kg here too)
seabirds$NitrogenInput.perArea <- (0.181 * (seabirds$Defecation.rate/1000) * seabirds$Estimate.adjusted * seabirds$Season.length)/
  seabirds$Area
seabirds$NitrogenInput.perArea.low <- (0.181 * (seabirds$Defecation.rate/1000) * seabirds$Estimate.adjusted.low * seabirds$Season.length)/
  seabirds$Area
seabirds$NitrogenInput.perArea.high <- (0.181 * (seabirds$Defecation.rate/1000) * seabirds$Estimate.adjusted.high * seabirds$Season.length)/
  seabirds$Area

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Total nitrogen input:
seabirds$tot.Nitrogen <- seabirds$NitrogenInput.perArea * seabirds$Area
seabirds$tot.Nitrogen.low <- seabirds$NitrogenInput.perArea.low * seabirds$Area
seabirds$tot.Nitrogen.high <- seabirds$NitrogenInput.perArea.high * seabirds$Area

write.csv(seabirds, "SNP_Data/Processed/Seabird_NutrientInput_Predicted_OtherAtolls.csv")

# In tonnes:
colSums(seabirds[,c("tot.Nitrogen", "tot.Nitrogen.low", "tot.Nitrogen.high")])/1000
