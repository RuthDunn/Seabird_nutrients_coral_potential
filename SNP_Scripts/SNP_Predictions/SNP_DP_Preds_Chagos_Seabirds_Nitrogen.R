rm(list = ls(all = TRUE))

library(plyr)
library(ggplot2)
library(ape)
library(data.table)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Data loads ####

# Load Chagos seabird census data:

seabirds <- read.csv("SNP_Data/Carr_etal_2020_Chagos_breeding_seabirds.csv")
names(seabirds)[1] <- "Island"
seabirds$Atoll_Island <- paste(seabirds$Atoll, seabirds$Island, sep = "_")
seabirds <- seabirds[,-c(1:3)]
seabirds <- reshape2::melt(data = seabirds, 
                           id.vars = c("Atoll_Island"),
                           variable.name = "Species",
                           value.name = "Census")

# Load predicted seabird numbers:
pred.high.nn <- read.csv("SNP_ModelOutputs/Chagos_Seabird_Predictions_HighNNForest_priors.csv")
pred.high.nn <- pred.high.nn[,c("Atoll_Island", "Species", "Estimate", "Q2.5", "Q97.5")]
pred.low.nn <- read.csv("SNP_ModelOutputs/Chagos_Seabird_Predictions_LowNNForest_priors.csv")
pred.low.nn <- pred.low.nn[,c("Atoll_Island", "Species", "Estimate", "Q2.5", "Q97.5")]
pred.mid.nn <- read.csv("SNP_ModelOutputs/Chagos_Seabird_Predictions_MidNNForest_priors.csv")
pred.mid.nn <- pred.mid.nn[,c("Atoll_Island", "Species", "Estimate", "Q2.5", "Q97.5")]

# Add predicted values
pred.high.nn <- join(seabirds, pred.high.nn, by = c("Atoll_Island", "Species"))
pred.high.nn$Value.high.nn <- ifelse(is.na(pred.high.nn$Estimate) == T, pred.high.nn$Census, pred.high.nn$Estimate)
pred.high.nn$Value.high.nn.q2.5 <- ifelse(is.na(pred.high.nn$Q2.5) == T, pred.high.nn$Census, pred.high.nn$Q2.5)
pred.high.nn$Value.high.nn.q97.5 <- ifelse(is.na(pred.high.nn$Q97.5) == T, pred.high.nn$Census, pred.high.nn$Q97.5)
pred.high.nn <- pred.high.nn[,c("Atoll_Island", "Species", "Value.high.nn", "Value.high.nn.q2.5", "Value.high.nn.q97.5")]

pred.low.nn <- join(seabirds, pred.low.nn, by = c("Atoll_Island", "Species"))
pred.low.nn$Value.low.nn <- ifelse(is.na(pred.low.nn$Estimate) == T, pred.low.nn$Census, pred.low.nn$Estimate)
pred.low.nn$Value.low.nn.q2.5 <- ifelse(is.na(pred.low.nn$Q2.5) == T, pred.low.nn$Census, pred.low.nn$Q2.5)
pred.low.nn$Value.low.nn.q.97.5 <- ifelse(is.na(pred.low.nn$Q97.5) == T, pred.low.nn$Census, pred.low.nn$Q97.5)
pred.low.nn <- pred.low.nn[,c("Atoll_Island", "Species", "Value.low.nn", "Value.low.nn.q2.5", "Value.low.nn.q.97.5")]

pred.mid.nn <- join(seabirds, pred.mid.nn, by = c("Atoll_Island", "Species"))
pred.mid.nn$Value.mid.nn <- ifelse(is.na(pred.mid.nn$Estimate) == T, pred.mid.nn$Census, pred.mid.nn$Estimate)
pred.mid.nn$Value.mid.nn.q2.5 <- ifelse(is.na(pred.mid.nn$Q2.5) == T, pred.mid.nn$Census, pred.mid.nn$Q2.5)
pred.mid.nn$Value.mid.nn.q.97.5 <- ifelse(is.na(pred.mid.nn$Q97.5) == T, pred.mid.nn$Census, pred.mid.nn$Q97.5)
pred.mid.nn <- pred.mid.nn[,c("Atoll_Island", "Species", "Value.mid.nn", "Value.mid.nn.q2.5", "Value.mid.nn.q.97.5")]

seabirds <- cbind(seabirds,
                  pred.high.nn[,c("Value.high.nn", "Value.high.nn.q2.5", "Value.high.nn.q97.5")],
                  pred.low.nn[,c("Value.low.nn", "Value.low.nn.q2.5", "Value.low.nn.q.97.5")],
                  pred.mid.nn[,c("Value.mid.nn", "Value.mid.nn.q2.5", "Value.mid.nn.q.97.5")])

rm(pred.high.nn, pred.low.nn, pred.mid.nn)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Nitrogen input per ha (NI) calculation:

# (Graham et al 2018)
# NI[i,j] = N[g] * Dr[i] * Bd[i,j] * Res[i,j] / IsArea[j]

# N[g] = nitrogen content of guano = 18.1% 
# Dr[i] = daily defecation rate (g) per species i ~ scaled allometrically from boobies
# Res = number of days of the year that the species is resident
# IsArea = island area

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Bring in:
#       - Island area
#       - Seabird breeding season length
#       - Seabird mass

# 1) Island area:
areas <- read.csv("SNP_Data/Carr_2021_Chagos_islands_info.csv")
areas$Atoll_Island <- as.factor(paste(areas$Atoll, areas$Island, sep = "_"))
# Revalue Peros_Banhos_Anglaise (from Peros_Banhos_Anglais):
levels(areas$Atoll_Island)[levels(areas$Atoll_Island) == "Peros_Banhos_Anglais"] <- "Peros_Banhos_Anglaise"
seabirds <- join(seabirds, areas[,c("Atoll_Island", "Size_Ha")], by = "Atoll_Island")
# seabirds <- seabirds[complete.cases(seabirds), ]
rm(areas)

# 2) Seabird breeding season length
# 3) Seabird mass
bird.info <- read.csv("SNP_Data/Schreiber_Burger_2002_Seabird_info.csv")
names(bird.info)[1] <- "Species"
bird.info$Season.length <- bird.info$IncPhase + bird.info$BroodPhase + bird.info$CrechePhase
seabirds <- join(seabirds, bird.info[,c("Species", "Season.length", "Mass")], by = "Species")
rm(bird.info)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Do some kind of guano calculation
# "We found birds defecated approximately every 36 min, each dropping weighing 0.66 g"
# So daily defecation rate = 
# (1440/36) * 0.66 = 26.4 g/day
# "For calculating defecation quantities for other species, we assume guano produced scaled
# allometrically with body size of species." (Young et al., 2010)
# DR[i] = DR[booby] * Mass[i]^0.63/Mass[booby]^0.63 (Staunton Smith & Johnson 1995)

seabirds$Defecation.rate <- (26.4 * (seabirds$Mass*1000)^0.63)/1150^0.63

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# NI[i,j] = N[g] * Dr[i] * Bd[i,j] * Res[i,j]
# Include conversion to kg/year:

# Adjust Bd (number of days on the island) to account for time away whilst foraging
# Sooty tern - one bird foraging - assume one adult present
# Red-footed booby - one adult always present, one away over night - assume 1.5 birds present
# Lesser noddy - short foraging trips, don't adjust numbers - assume 2 birds present

seabirds <- seabirds %>%
  mutate(Census = ifelse(Species == "Sula_sula", Census * 1.5,
                         ifelse(Species == "Anous_tenuirostris", Census * 2,
                                Census)),
         Value.high.nn = ifelse(Species == "Sula_sula", Value.high.nn * 1.5,
                         ifelse(Species == "Anous_tenuirostris", Value.high.nn * 2,
                                Value.high.nn)),
         Value.low.nn = ifelse(Species == "Sula_sula", Value.low.nn * 1.5,
                         ifelse(Species == "Anous_tenuirostris", Value.low.nn * 2,
                                Value.low.nn)))

seabirds$NitrogenInput.Current <- (0.181 * (seabirds$Defecation.rate/1000) * seabirds$Census * seabirds$Season.length)
seabirds$NitrogenInput.high.nn <- (0.181 * (seabirds$Defecation.rate/1000) * seabirds$Value.high.nn * seabirds$Season.length)
seabirds$NitrogenInput.low.nn <- (0.181 * (seabirds$Defecation.rate/1000) * seabirds$Value.low.nn * seabirds$Season.length)
seabirds$NitrogenInput.mid.nn <- (0.181 * (seabirds$Defecation.rate/1000) * seabirds$Value.mid.nn * seabirds$Season.length)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Check that this all looks okay:

library(tidyverse)

test <- seabirds %>%
  group_by(Atoll_Island) %>%
  summarise(CurrentN = sum(NitrogenInput.Current),
            Lownn.N = sum(NitrogenInput.low.nn),
            Highnn.N = sum(NitrogenInput.high.nn),
            Midnn.N = sum(NitrogenInput.mid.nn),
            Area = mean(Size_Ha))

sum(test$CurrentN)/1000
sum(test$Highnn.N)/1000
sum(test$Lownn.N)/1000
sum(test$Midnn.N)/1000

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# We're going to use this dataframe to make predictions
# Therefore, keep "rat islands" only (other than DG):

rat.islands <- read.csv("SNP_Data/Processed/ChagosRatIslands.csv")
rat.islands <- rat.islands[-1,c("Atoll_Island")]

# Select these rows only:
seabirds.rats <- seabirds[is.element(seabirds$Atoll_Island, rat.islands),]

nlevels(as.factor(seabirds.rats$Atoll_Island))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

write.csv(seabirds.rats, "SNP_Data/Processed/Seabird_NutrientInput_Chagos_Predicted_RatIslandsOnly_priors.csv")
write.csv(seabirds, "SNP_Data/Processed/Seabird_NutrientInput_Chagos_Predicted_priors.csv")
