rm(list = ls(all = TRUE))

# Packages:
library(plyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(ggradar)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Prepare seabird population prediction data ####

# 1) Current data

current.birds <- read.csv("SNP_Data/Processed/SeabirdDensityDrivers.csv")
current.birds <- current.birds[,c("Atoll_Island", "Rattus_rattus", "Size_Ha", "NonNativeForest_p",
                                  "Sula_sula", "Onychoprion_fuscatus", "Anous_tenuirostris")]
# Ratty islands only:
current.birds <- current.birds[which(current.birds$Rattus_rattus=='P'),]
# Islands with area only:
current.birds <- current.birds[!is.na(current.birds$Size_Ha),]
# Remove DG main island:
current.birds <- current.birds[which(current.birds$Atoll_Island != "Diego_Garcia_Diego_Garcia"),]
# Create vector of atoll_islands:
atoll_islands <- unique(current.birds$Atoll_Island)

# Summarise:
Estimate <- colSums(current.birds[,c(5:7)], dims = 1)
Species <- colnames(current.birds[,5:7])
Current <- as.data.frame(cbind(Species, Estimate))
rm(current.birds, Estimate, Species)
Current$High <- NA
Current$Low <- NA
Current$Type <- "Current"

# ~~~~~~~~~~~~~~~
# 2) Prediction data
# a) High non-native forest cover

# Load data:
pred.high.nn <- read.csv("SNP_ModelOutputs/Chagos_Seabird_Predictions_HighNNForest.hln.csv")
# Select cols
pred.high.nn <- pred.high.nn[,c("Atoll_Island", "Species", "Estimate", "Q2.5", "Q97.5")]
# Rename cols
colnames(pred.high.nn) <- c("Atoll_Island", "Species", "Estimate.high.nn", "Q2.5.high.nn", "Q97.5.high.nn")
# Transform logged data
# pred.high.nn$Estimate.high.nn <- exp(pred.high.nn$Estimate.high.nn)
# pred.high.nn$Q2.5.high.nn <- exp(pred.high.nn$Q2.5.high.nn)
# pred.high.nn$Q97.5.high.nn <- exp(pred.high.nn$Q97.5.high.nn)
# Keep islands that we have current bird estimates for only:
pred.high.nn <- pred.high.nn[which(pred.high.nn$Atoll_Island == "Great_Chagos_Bank_Eagle_Island"|
                                     pred.high.nn$Atoll_Island == "Peros_Banhos_Coin"|
                                     pred.high.nn$Atoll_Island == "Peros_Banhos_Monpatre_complex"|
                                     pred.high.nn$Atoll_Island == "Peros_Banhos_Poule"|
                                     pred.high.nn$Atoll_Island == "Peros_Banhos_Petit_Souer"|
                                     pred.high.nn$Atoll_Island == "Peros_Banhos_Grand_Souer"|
                                     pred.high.nn$Atoll_Island == "Peros_Banhos_Pierre"|
                                     pred.high.nn$Atoll_Island == "Peros_Banhos_Petite_Mapou"|
                                     pred.high.nn$Atoll_Island == "Peros_Banhos_Grande_Mapou"|
                                     pred.high.nn$Atoll_Island == "Peros_Banhos_Diamant"|
                                     pred.high.nn$Atoll_Island == "Peros_Banhos_Passe"|
                                     pred.high.nn$Atoll_Island == "Peros_Banhos_Moresby"|
                                     pred.high.nn$Atoll_Island == "Peros_Banhos_Yeye"|
                                     pred.high.nn$Atoll_Island == "Peros_Banhos_Fouquet"|
                                     pred.high.nn$Atoll_Island == "Peros_Banhos_Mapou_de_Coin"|
                                     pred.high.nn$Atoll_Island == "Salomon_Islands_Boddam"|
                                     pred.high.nn$Atoll_Island == "Salomon_Islands_Diable"|
                                     pred.high.nn$Atoll_Island == "Salomon_Islands_Anglaise"|
                                     pred.high.nn$Atoll_Island == "Salomon_Islands_Takamaka"|
                                     pred.high.nn$Atoll_Island == "Salomon_Islands_Fouquet"|
                                     pred.high.nn$Atoll_Island == "Salomon_Islands_Sepulture"|
                                     pred.high.nn$Atoll_Island == "Salomon_Islands_Poule"),]

# Summarise:
High.nn.est <- aggregate(Estimate.high.nn ~ Species, pred.high.nn, sum)
High.nn.high <- aggregate(Q97.5.high.nn ~ Species, pred.high.nn, sum)
High.nn.low <- aggregate(Q2.5.high.nn ~ Species, pred.high.nn, sum)
High.nn <- cbind(High.nn.est, High.nn.high[,2], High.nn.low[,2])
colnames(High.nn) <- c("Species", "Estimate", "High", "Low")
High.nn$Type <- "HighNN"
rm(pred.high.nn, High.nn.est, High.nn.low, High.nn.high)

# ~~~~~~~~~~~~~~~
# b) Low non-native forest cover

pred.low.nn <- read.csv("SNP_ModelOutputs/Chagos_Seabird_Predictions_LowNNForest.hln.csv")
pred.low.nn <- pred.low.nn[,c("Atoll_Island", "Species", "Estimate", "Q2.5", "Q97.5")]
colnames(pred.low.nn) <- c("Atoll_Island", "Species", "Estimate.low.nn", "Q2.5.low.nn", "Q97.5.low.nn")
pred.low.nn$Atoll_Island_Sp <- paste(pred.low.nn$Atoll_Island, pred.low.nn$Species, sep = "_")
# Keep islands that we have current bird estimates for only:
pred.low.nn <- pred.low.nn[which(pred.low.nn$Atoll_Island == "Great_Chagos_Bank_Eagle_Island"|
                                   pred.low.nn$Atoll_Island == "Peros_Banhos_Coin"|
                                   pred.low.nn$Atoll_Island == "Peros_Banhos_Monpatre_complex"|
                                   pred.low.nn$Atoll_Island == "Peros_Banhos_Poule"|
                                   pred.low.nn$Atoll_Island == "Peros_Banhos_Petit_Souer"|
                                   pred.low.nn$Atoll_Island == "Peros_Banhos_Grand_Souer"|
                                   pred.low.nn$Atoll_Island == "Peros_Banhos_Pierre"|
                                   pred.low.nn$Atoll_Island == "Peros_Banhos_Petite_Mapou"|
                                   pred.low.nn$Atoll_Island == "Peros_Banhos_Grande_Mapou"|
                                   pred.low.nn$Atoll_Island == "Peros_Banhos_Diamant"|
                                   pred.low.nn$Atoll_Island == "Peros_Banhos_Passe"|
                                   pred.low.nn$Atoll_Island == "Peros_Banhos_Moresby"|
                                   pred.low.nn$Atoll_Island == "Peros_Banhos_Yeye"|
                                   pred.low.nn$Atoll_Island == "Peros_Banhos_Fouquet"|
                                   pred.low.nn$Atoll_Island == "Peros_Banhos_Mapou_de_Coin"|
                                   pred.low.nn$Atoll_Island == "Salomon_Islands_Boddam"|
                                   pred.low.nn$Atoll_Island == "Salomon_Islands_Diable"|
                                   pred.low.nn$Atoll_Island == "Salomon_Islands_Anglaise"|
                                   pred.low.nn$Atoll_Island == "Salomon_Islands_Takamaka"|
                                   pred.low.nn$Atoll_Island == "Salomon_Islands_Fouquet"|
                                   pred.low.nn$Atoll_Island == "Salomon_Islands_Sepulture"|
                                   pred.low.nn$Atoll_Island == "Salomon_Islands_Poule"),]

# Summarise:
Est <- aggregate(Estimate.low.nn ~ Species, pred.low.nn, sum)
High <- aggregate(Q97.5.low.nn ~ Species, pred.low.nn, sum)
Low <- aggregate(Q2.5.low.nn ~ Species, pred.low.nn, sum)
Low.nn <- cbind(Est, High[,2], Low[,2])
colnames(Low.nn) <- c("Species", "Estimate", "High", "Low")
Low.nn$Type <- "LowNN"
rm(pred.low.nn, Est, High, Low)

seabird.pops <- rbind(Current, High.nn, Low.nn)
rm(High.nn, Low.nn, Current, atoll_islands)

seabird.pops$Estimate <- as.numeric(seabird.pops$Estimate)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Prepare prey available/seabird consumption data ####

# 1) Load consumption data

consumption <- read.csv("SNP_Data/Processed/Seabird_Group_Consumption_Predicted.csv")
names(consumption)[1] <- "Species"

# Reorder:
consumption <- consumption[,c("Species",
                              "Pop_Consumption_t", "Pop_Consumption_t.pred.low.nn", "Pop_Consumption_t.pred.high.nn",
                              "Pop_Consumption_low_t", "Pop_Consumption_t.pred.low.nn.low", "Pop_Consumption_t.pred.high.nn.low",
                              "Pop_Consumption_high_t", "Pop_Consumption_t.pred.low.nn.high", "Pop_Consumption_t.pred.high.nn.high"),]
# Wide to long:
consumption.vals <- gather(consumption[,c(1:4)], Biomass, Estimate, Pop_Consumption_t:Pop_Consumption_t.pred.high.nn)
consumption.low <- gather(consumption[,c(1, 5:7)], Biomass, Low, Pop_Consumption_low_t:Pop_Consumption_t.pred.high.nn.low)
consumption.high <- gather(consumption[,c(1, 8:10)], Biomass, High, Pop_Consumption_high_t:Pop_Consumption_t.pred.high.nn.high)

consumption.vals$Low <- consumption.low$Low
consumption.vals$High <- consumption.high$High
consumption <- consumption.vals

rm(consumption.vals, consumption.low, consumption.high)

consumption$dx <- rep(c(0.2, -0.2, 0), each = 4)
consumption$x0 <- c(1,2,3,4)

# 2) Load prey data

prey <- read.csv("SNP_Data/Processed/OceanBiomass_BufferAreas_R.csv")
prey$X <- c("Sula_sula", "Anous_tenuirostris", "Onychoprion_fuscatus", "All")
names(prey)[1] <- "Species"
prey <- prey[,c("Species", "meantonnes", "sdtonnes")]
prey$Biomass <- "Available"
names(prey)[2] <- "Estimate"
prey$Lower <- prey$Estimate-prey$sdtonnes
prey$Upper <- prey$Estimate+prey$sdtonnes
prey$dx <- -0.3
prey$x0 <- c(4,2,3,1)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Plot seabird population predictions ####

p1 <- ggplot(seabird.pops, aes(x = Species, y = (Estimate+1), ymin = (Low+1), ymax = (High+1), fill = Type)) +
  geom_bar(stat = "identity",
           width = .6, position = position_dodge()) +
  scale_fill_manual(values = c("#DC267F", "#785EF0", "#648FFF"),
                    labels = c("Current", "Rat eradication", "Rat eradication & habitat restoration")) +
  
  geom_errorbar(stat = "identity",
                width = .3, position = position_dodge(0.6), aes(col = Type)) +
  scale_colour_manual(values = c("#DC267F", "#4e34c7", "#3d61bf")) +
  
  theme_light() %+replace% theme(panel.grid.minor = element_blank(),
                                 legend.position = 'bottom') +
  ylab("Seabird abundance (number of pairs)") +
  scale_x_discrete(labels = c("Lesser noddy", "Sooty tern", "Red-footed booby")) +
  guides(fill = guide_legend(title = element_text("Scenario")), col = "none") +
  scale_y_continuous(labels = scales::comma, trans = "log10")

p1

# ggsave("Plots/Seabird_Predictions_ppt.png", width = 5, height = 4)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Plot seabird consumption predictions ####

p2 <- ggplot() +
  geom_bar(data = prey,
           aes(y = rev(Species), x = Estimate), stat = "identity",
           fill = "#FFB000", alpha = .4, width = 1) +
  geom_bar(data = consumption,
           aes(y = (x0 + dx), x = (Estimate), fill = Biomass),
           stat = "identity", width = .2) +
  geom_errorbar(data = consumption, aes(y = (x0 + dx), xmin = Low, xmax = High, col = Biomass),
                stat = "identity", width = .1) +
  scale_fill_manual(values = c("#DC267F", "#785EF0", "#648FFF")) +
  # labels = c("Current", "Rat eradication", "Rat eradication & habitat restoration")) +
  scale_colour_manual(values = c("#ed7bb3", "#4e34c7", "#3d61bf")) +
  ylab(element_blank()) +
  xlab("Seabird prey requirements (t)") +
  theme_light() %+replace% theme(panel.grid.minor = element_blank(),
                                 legend.position = 'bottom', legend.title = element_blank()) +
  guides(col = "none") +
  scale_x_continuous(labels = scales::comma, trans = "log10") +
  scale_y_discrete(labels = c("Red-footed booby", "Sooty tern", "Lesser noddy", "All species")) +
  annotate("text", x = 800000, y = 4, label = "Available biomass", angle = 270, col = "#FFB000", size = 3)

p2

p3 <- ggplot() +
  geom_bar(data = prey,
           aes(x = rev(Species), y = Estimate), stat = "identity",
           fill = "#FFB000", alpha = .4, width = 1) +
  geom_bar(data = consumption,
           aes(x = (x0 + dx), y = (Estimate), fill = Biomass),
           stat = "identity", width = .2) +
  geom_errorbar(data = consumption, aes(x = (x0 + dx), ymin = Low, ymax = High, col = Biomass),
                stat = "identity", width = .1) +
  scale_fill_manual(values = c("#DC267F", "#785EF0", "#648FFF")) +
  # labels = c("Current", "Rat eradication", "Rat eradication & habitat restoration")) +
  scale_colour_manual(values = c("#ed7bb3", "#4e34c7", "#3d61bf")) +
  xlab(element_blank()) +
  ylab("Seabird prey requirements (t)") +
  theme_light() %+replace% theme(panel.grid.minor = element_blank(),
                                 legend.position = 'bottom', legend.title = element_blank()) +
  guides(col = "none") +
  scale_y_continuous(labels = scales::comma, trans = "log10") +
  scale_x_discrete(labels = c("Red-footed booby", "Sooty tern", "Lesser noddy", "All species")) +
  annotate("text", y = 8000, x = 1, label = "Available biomass", col = "#FFB000", size = 3)

p3

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Save plot for manuscript ####

# For manuscript:

ggarrange(p1, p2,  common.legend = T, legend = "bottom", widths = c(1.5,1))

ggsave("Plots/Seabird_Predictions_ms.png", width = 8, height = 5)

# For poster

ggarrange(p1, p3,  common.legend = T, legend = "bottom")

ggsave("Plots/Seabird_Predictions_poster.png", width = 55, height = 15, units = "cm")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Numbers for manuscripts ####

# 1) Seabird pop changes

# Current pops:
seabird.pops[c("Onychoprion_fuscatus", "Anous_tenuirostris", "Sula_sula"), "Estimate"]

# Pop change with eradication:
seabird.pops[c("1", "2", "3"),]

# Pop change ft eradication & habitat restoration:
# LN:
((seabird.pops["21", "Estimate"] - seabird.pops["2","Estimate"])/
    seabird.pops["2","Estimate"]) * 100
# ST:
((seabird.pops["11", "Estimate"] - seabird.pops["1","Estimate"])/
    seabird.pops["1","Estimate"]) * 100
# RFB:
((seabird.pops["31", "Estimate"] - seabird.pops["3","Estimate"])/
    seabird.pops["3","Estimate"]) * 100

# Total % change:
((seabird.pops["21", "Estimate"] + seabird.pops["11", "Estimate"] + seabird.pops["31", "Estimate"]) -
  (seabird.pops["2","Estimate"] + seabird.pops["1","Estimate"] + seabird.pops["3","Estimate"]))/
  (seabird.pops["2","Estimate"] + seabird.pops["1","Estimate"] + seabird.pops["3","Estimate"]) *
  100

seabird.pops["21", "Estimate"] + seabird.pops["11", "Estimate"] + seabird.pops["31", "Estimate"]
seabird.pops["21", "Low"] + seabird.pops["11", "Low"] + seabird.pops["31", "Low"]
seabird.pops["21", "High"] + seabird.pops["11", "High"] + seabird.pops["31", "High"]
