rm(list = ls(all = TRUE))

library(plyr)
library(ggplot2)
library(ape)
library(MCMCglmm)
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

# Add predicted values
pred.high.nn <- join(seabirds, pred.high.nn, by = c("Atoll_Island", "Species"))
pred.high.nn$Value.high.nn <- ifelse(is.na(pred.high.nn$Estimate) == T, pred.high.nn$Census, pred.high.nn$Estimate)
pred.high.nn$Value.high.nn.low <- ifelse(is.na(pred.high.nn$Estimate) == T, pred.high.nn$Census, pred.high.nn$Q2.5)
pred.high.nn$Value.high.nn.high <- ifelse(is.na(pred.high.nn$Estimate) == T, pred.high.nn$Census, pred.high.nn$Q97.5)

pred.low.nn <- join(seabirds, pred.low.nn, by = c("Atoll_Island", "Species"))
pred.low.nn$Value.low.nn <- ifelse(is.na(pred.low.nn$Estimate) == T, pred.low.nn$Census, pred.low.nn$Estimate)
pred.low.nn$Value.low.nn.low <- ifelse(is.na(pred.low.nn$Estimate) == T, pred.low.nn$Census, pred.low.nn$Q2.5)
pred.low.nn$Value.low.nn.high <- ifelse(is.na(pred.low.nn$Estimate) == T, pred.low.nn$Census, pred.low.nn$Q97.5)

seabirds <- cbind(seabirds,
                  pred.high.nn[,c("Value.high.nn", "Value.high.nn.low", "Value.high.nn.high")],
                  pred.low.nn[,c("Value.low.nn", "Value.low.nn.low", "Value.low.nn.high")])

rm(pred.high.nn, pred.low.nn)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Load seabird masses and breeding stage lengths ####
# (from Biology of Marine Birds, Schreiber and Burger, 2002, Appendix 2

seabirds.info <- read.csv("SNP_Data/Schreiber_Burger_2002_Seabird_info.csv")
# Change mass from kg to g
seabirds.info$Mass <- seabirds.info$Mass * 1000
names(seabirds.info)[1] <- "Species"

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Load Chagos island latitude data ####

areas <- read.csv("SNP_Data/Carr_2021_Chagos_islands_info.csv")
names(areas)[1] <- "Island"
areas$Atoll_Island <- paste(areas$Atoll, areas$Island, sep = "_")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Join dfs together ####

seabirds <- join(seabirds, seabirds.info[,c(1:5)], by = "Species")
seabirds <- join(seabirds, areas, by = "Atoll_Island")
rm(areas, seabirds.info)

# Remove NA rows:
seabirds<-seabirds[complete.cases(seabirds),]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Load in seabird FMR calculator model :) ####

load(file = "C:/Users/Ruth/Dropbox/PhD/Project/Breeding Energetics - Meta Analysis/FMR_meta_analysis/myapp/data/model4.rda")

# Create new column of bird species that match those in the model
seabirds$Species_Model <- paste("animal.", seabirds$Species, sep = "")

seabirds$Species_Model <- revalue(seabirds$Species_Model, c("animal.Ardenna_pacifica" = "animal.Puffinus_pacificus"))
seabirds$Species_Model <- revalue(seabirds$Species_Model, c("animal.Fregeta_ariel" = "animal.Fregata_ariel"))
seabirds$Species_Model <- revalue(seabirds$Species_Model, c("animal.Puffinus_bailloni" = "animal.Puffinus_lherminieri"))
seabirds$Species_Model <- revalue(seabirds$Species_Model, c("animal.Fregeta_minor" = "animal.Fregata_minor"))
seabirds$Species_Model <- revalue(seabirds$Species_Model, c("animal.Thalasseus_bergii" = "animal.Sterna_bergii"))
seabirds$Species_Model <- revalue(seabirds$Species_Model, c("animal.Sternula_albifrons" = "animal.Sterna_albifrons"))
seabirds$Species_Model <- revalue(seabirds$Species_Model, c("animal.Onychoprion_anaethetus" = "animal.Sterna_anaethetus"))
seabirds$Species_Model <- revalue(seabirds$Species_Model, c("animal.Onychoprion_fuscatus" = "animal.Sterna_fuscata"))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Make estimates of FMR ####

# 1) Incubation

# Incubation posterior
NewVal <- NULL
for(i in 1:nrow(seabirds)){
  # i = 5
  NewVal_i <- 10^posterior.mode(model4$Sol[,"(Intercept)"] +
                                           model4$Sol[,seabirds$Species_Model[i]] +
                                           seabirds$Centroid_Lat[i] * model4$Sol[,"Lat"] +
                                           log10(seabirds$Mass[i])*model4$Sol[,"log_Mass"] +
                                           model4$Sol[,"PhaseIncubation"])
  NewVal <- rbind(NewVal_i, NewVal)
}
seabirds$FMR_Inc_posterior <- as.numeric(NewVal)

# Incubation HPD interval
NewVal <- NULL
for(i in 1:nrow(seabirds)){
  # i = 278
  NewVal_i <- 10^HPDinterval(model4$Sol[,"(Intercept)"] +
                                          model4$Sol[,seabirds$Species_Model[i]] +
                                          seabirds$Centroid_Lat[i] * model4$Sol[,"Lat"] +
                                          log10(seabirds$Mass[i])*model4$Sol[,"log_Mass"] +
                                          model4$Sol[,"PhaseIncubation"])
  NewVal <- rbind(NewVal_i, NewVal)
}
seabirds$FMR_Inc_HPD_lower <- NewVal[,1]
seabirds$FMR_Inc_HPD_upper <- NewVal[,2]

# 2) Brood

# Brood posterior
NewVal <- NULL
for(i in 1:nrow(seabirds)){
  # i = 278
  NewVal_i <- 10^posterior.mode(model4$Sol[,"(Intercept)"] +
                                  model4$Sol[,seabirds$Species_Model[i]] +
                                  seabirds$Centroid_Lat[i] * model4$Sol[,"Lat"] +
                                  log10(seabirds$Mass[i])*model4$Sol[,"log_Mass"])
  NewVal <- rbind(NewVal_i, NewVal)
}
seabirds$FMR_Brood_posterior <- as.numeric(NewVal)

# Brood HPD interval
NewVal <- NULL
for(i in 1:nrow(seabirds)){
  # i = 278
  NewVal_i <- 10^HPDinterval(model4$Sol[,"(Intercept)"] +
                               model4$Sol[,seabirds$Species_Model[i]] +
                               seabirds$Centroid_Lat[i] * model4$Sol[,"Lat"] +
                               log10(seabirds$Mass[i])*model4$Sol[,"log_Mass"])
  NewVal <- rbind(NewVal_i, NewVal)
}
seabirds$FMR_Brood_HPD_lower <- NewVal[,1]
seabirds$FMR_Brood_HPD_upper <- NewVal[,2]

# 3) Creche

# Creche posterior
NewVal <- NULL
for(i in 1:nrow(seabirds)){
  # i = 278
  NewVal_i <- 10^posterior.mode(model4$Sol[,"(Intercept)"] +
                                  model4$Sol[,seabirds$Species_Model[i]] +
                                  seabirds$Centroid_Lat[i] * model4$Sol[,"Lat"] +
                                  log10(seabirds$Mass[i])*model4$Sol[,"log_Mass"] +
                                  model4$Sol[,"PhaseCreche"])
  NewVal <- rbind(NewVal_i, NewVal)
}
seabirds$FMR_Creche_posterior <- as.numeric(NewVal)

# Creche HPD interval
NewVal <- NULL
for(i in 1:nrow(seabirds)){
  # i = 278
  NewVal_i <- 10^HPDinterval(model4$Sol[,"(Intercept)"] +
                               model4$Sol[,seabirds$Species_Model[i]] +
                               seabirds$Centroid_Lat[i] * model4$Sol[,"Lat"] +
                               log10(seabirds$Mass[i])*model4$Sol[,"log_Mass"] +
                               model4$Sol[,"PhaseIncubation"])
  NewVal <- rbind(NewVal_i, NewVal)
}
seabirds$FMR_Creche_HPD_lower <- NewVal[,1]
seabirds$FMR_Creche_HPD_upper <- NewVal[,2]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Calc total breeding season FMR ####

seabirds$Total_FMR <- seabirds$IncPhase * seabirds$FMR_Inc_posterior +
  seabirds$BroodPhase * seabirds$FMR_Brood_posterior +
  seabirds$CrechePhase * seabirds$FMR_Creche_posterior

seabirds$Total_FMR_low <- seabirds$IncPhase * seabirds$FMR_Inc_HPD_lower +
  seabirds$BroodPhase * seabirds$FMR_Brood_HPD_lower +
  seabirds$CrechePhase * seabirds$FMR_Creche_HPD_lower

seabirds$Total_FMR_up <- seabirds$IncPhase * seabirds$FMR_Inc_HPD_upper +
  seabirds$BroodPhase * seabirds$FMR_Brood_HPD_upper +
  seabirds$CrechePhase * seabirds$FMR_Creche_HPD_upper

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Convert FMR to biomass consumption ####
# (Take into account assimilation rate and energy content of prey)

seabirds$Indiv_Consumption_g <- (0.75 * seabirds$Total_FMR)/5.5
seabirds$Indiv_Consumption_low_g <- (0.75 * seabirds$Total_FMR_low)/5.5
seabirds$Indiv_Consumption_up_g <- (0.75 * seabirds$Total_FMR_up)/5.5

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Multiply by population size ####
# (which is currently in pairs, so multiply by 2)
# Also convert from grams into tonnes

seabirds$Pop_Consumption_t <- seabirds$Indiv_Consumption_g * (seabirds$Census * 2)/1000000
seabirds$Pop_Consumption_low_t <- seabirds$Indiv_Consumption_low_g * (seabirds$Census * 2)/1000000
seabirds$Pop_Consumption_high_t <- seabirds$Indiv_Consumption_up_g * (seabirds$Census * 2)/1000000

seabirds$Pop_Consumption_t.pred.low.nn <- seabirds$Indiv_Consumption_g * (seabirds$Value.low.nn * 2)/1000000
seabirds$Pop_Consumption_t.pred.low.nn.low <- seabirds$Indiv_Consumption_low_g * (seabirds$Value.low.nn * 2)/1000000
seabirds$Pop_Consumption_t.pred.low.nn.high <- seabirds$Indiv_Consumption_up_g * (seabirds$Value.low.nn * 2)/1000000

seabirds$Pop_Consumption_t.pred.high.nn <- seabirds$Indiv_Consumption_g * (seabirds$Value.high.nn * 2)/1000000
seabirds$Pop_Consumption_t.pred.high.nn.low <- seabirds$Indiv_Consumption_low_g * (seabirds$Value.high.nn * 2)/1000000
seabirds$Pop_Consumption_t.pred.high.nn.high <- seabirds$Indiv_Consumption_up_g * (seabirds$Value.high.nn * 2)/1000000

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Group by species
Pop_Consumption_t <- by(seabirds$Pop_Consumption_t, list(seabirds$Species), FUN = sum)
Pop_Consumption_low_t <- by(seabirds$Pop_Consumption_low_t, list(seabirds$Species), FUN = sum)
Pop_Consumption_high_t <- by(seabirds$Pop_Consumption_high_t, list(seabirds$Species), FUN = sum)

Pop_Consumption_t.pred.low.nn <- by(seabirds$Pop_Consumption_t.pred.low.nn, list(seabirds$Species), FUN = sum)
Pop_Consumption_t.pred.low.nn.low <- by(seabirds$Pop_Consumption_t.pred.low.nn.low, list(seabirds$Species), FUN = sum)
Pop_Consumption_t.pred.low.nn.high <- by(seabirds$Pop_Consumption_t.pred.low.nn.high, list(seabirds$Species), FUN = sum)

Pop_Consumption_t.pred.high.nn <- by(seabirds$Pop_Consumption_t.pred.high.nn, list(seabirds$Species), FUN = sum)
Pop_Consumption_t.pred.high.nn.low <- by(seabirds$Pop_Consumption_t.pred.high.nn.low, list(seabirds$Species), FUN = sum)
Pop_Consumption_t.pred.high.nn.high <- by(seabirds$Pop_Consumption_t.pred.high.nn.high, list(seabirds$Species), FUN = sum)

sp.consumption <- as.data.frame(cbind(t(rbind(Pop_Consumption_t, Pop_Consumption_low_t, Pop_Consumption_high_t,
                                              Pop_Consumption_t.pred.low.nn, Pop_Consumption_t.pred.low.nn.low, Pop_Consumption_t.pred.low.nn.high,
                                              Pop_Consumption_t.pred.high.nn, Pop_Consumption_t.pred.high.nn.low, Pop_Consumption_t.pred.high.nn.high))))

sp.consumption[,1:9] <- sapply(sp.consumption[,1:9],as.numeric)

rm(Pop_Consumption_t, Pop_Consumption_low_t, Pop_Consumption_high_t,
   Pop_Consumption_t.pred.low.nn, Pop_Consumption_t.pred.low.nn.low, Pop_Consumption_t.pred.low.nn.high,
   Pop_Consumption_t.pred.high.nn, Pop_Consumption_t.pred.high.nn.low, Pop_Consumption_t.pred.high.nn.high,
   model4, NewVal, NewVal_i, i)

# Create "All species" column

All <- colSums(sp.consumption)
sp.consumption <- rbind(sp.consumption, All)
rownames(sp.consumption)[19] <- "All"

# Select key species:

sp.consumption <- sp.consumption[c("Sula_sula", "Onychoprion_fuscatus", "Anous_tenuirostris", "All"),]

# Save!

write.csv(sp.consumption, "SNP_Data/Processed/Seabird_Group_Consumption_Predicted.csv")

# Results text:
sp.consumption["Anous_tenuirostris",c(7:9)]
sp.consumption["Anous_tenuirostris",c(4:6)]

