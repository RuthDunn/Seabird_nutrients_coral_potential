rm(list = ls(all = TRUE))

library(plyr)
library(ggplot2)
library(ape)
library(MCMCglmm)
library(data.table)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Load Chagos seabird census data ####

seabirds <- read.csv("SNP_Data/Carr_etal_2020_Chagos_breeding_seabirds.csv")
names(seabirds)[1] <- "Island"
seabirds <- seabirds[,-c(2,3)]

colSums(seabirds[,c("Onychoprion_fuscatus", "Anous_tenuirostris", "Sula_sula")],)

# Make seabird df into long format:

seabirds <- reshape2::melt(data = seabirds, 
                             id.vars = c("Island"),
                             variable.name = "Species",
                             value.name = "Census")
# Remove 0s by first changing them to NAs
seabirds[seabirds==0] <- NA
seabirds<-seabirds[complete.cases(seabirds),]

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

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Join dfs together ####

seabirds <- join(seabirds, seabirds.info[,c(1:5)], by = "Species")
seabirds <- join(seabirds, areas, by = "Island")
rm(areas, seabirds.info)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Load in seabird FMR calculator model :) ####

load(file = "C:/Users/ruth-/Dropbox/PhD/Project/Breeding Energetics - Meta Analysis/FMR_meta_analysis/myapp/data/model4.rda")

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
  # i = 278
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
seabirds$Pop_Consumption_up_t <- seabirds$Indiv_Consumption_up_g * (seabirds$Census * 2)/1000000

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Group by species
Pop_Consumption_t <- by(seabirds$Pop_Consumption_t, list(seabirds$Species), FUN = sum)
Pop_Consumption_low_t <- by(seabirds$Pop_Consumption_low_t, list(seabirds$Species), FUN = sum)
Pop_Consumption_up_t <- by(seabirds$Pop_Consumption_up_t, list(seabirds$Species), FUN = sum)

Group <- c("Shearwater", "Shearwater", "Tropic/Booby/Frigate",  "Tropic/Booby/Frigate",
           "Tropic/Booby/Frigate",  "Tropic/Booby/Frigate",  "Tropic/Booby/Frigate",
           "Tropic/Booby/Frigate",  "Tropic/Booby/Frigate", "Tern/Noddy",
           "Tern/Noddy", "Tern/Noddy", "Tern/Noddy", "Tern/Noddy", "Tern/Noddy",
           "Tern/Noddy", "Tern/Noddy", "Tern/Noddy")

sp.consumption <- as.data.frame(cbind(t(rbind(Pop_Consumption_t, Pop_Consumption_low_t, Pop_Consumption_up_t)), Group))
sp.consumption[,1:3] <- sapply(sp.consumption[,1:3],as.numeric)


# Group by group

Pop_Consumption_t <- by(sp.consumption$Pop_Consumption_t, list(sp.consumption$Group), FUN = sum)
Pop_Consumption_low_t <- by(sp.consumption$Pop_Consumption_low_t, list(sp.consumption$Group), FUN = sum)
Pop_Consumption_up_t <- by(sp.consumption$Pop_Consumption_up_t, list(sp.consumption$Group), FUN = sum)

group.consumption <- as.data.frame(t(rbind(Pop_Consumption_t, Pop_Consumption_low_t, Pop_Consumption_up_t)))

# Add in total pop consumption

Pop_Consumption_t <- colSums(sp.consumption[1])
Pop_Consumption_low_t <- colSums(sp.consumption[2])
Pop_Consumption_up_t <- colSums(sp.consumption[3])

All <- cbind(Pop_Consumption_t, Pop_Consumption_low_t, Pop_Consumption_up_t)

group.consumption.all <- rbind(group.consumption, All)

write.csv(group.consumption.all, "SNP_Data/Processed/Seabird_Group_Consumption.csv")
