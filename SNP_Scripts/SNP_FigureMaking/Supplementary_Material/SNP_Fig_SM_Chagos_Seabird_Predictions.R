rm(list = ls(all = TRUE))

library(ggplot2)
library(RColorBrewer)
library(patchwork)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Load data ####

birds <- read.csv("SNP_Data/Processed/Seabird_NutrientInput_Chagos_Predicted.csv")
head(birds)
birds <- birds[,c(2:10)]

# Wide to long:

birds <- gather(birds, Estimate, Value, Census:Value.low.nn.q.97.5)

# Divide up data:

current <- birds[which(birds$Estimate == "Census"),]
high.nn <- birds[which(birds$Estimate == "Value.high.nn"),]
high.nn.q97.5 <- birds[which(birds$Estimate == "Value.high.nn.q97.5"),]
high.nn.q2.5 <- birds[which(birds$Estimate == "Value.high.nn.q2.5"),]
low.nn <- birds[which(birds$Estimate == "Value.low.nn"),]
low.nn.q97.5  <- birds[which(birds$Estimate == "Value.low.nn.q.97.5"),]
low.nn.q2.5 <- birds[which(birds$Estimate == "Value.low.nn.q2.5"),]

# Calculate all species sums:

current <- spread(current, Species, Value)
high.nn.q97.5 <- spread(high.nn.q97.5, Species, Value)
high.nn.q2.5 <- spread(high.nn.q2.5, Species, Value)
high.nn <- spread(high.nn, Species, Value)
low.nn <- spread(low.nn, Species, Value)
low.nn.q97.5 <- spread(low.nn.q97.5, Species, Value)
low.nn.q2.5 <- spread(low.nn.q2.5, Species, Value)

current$Current <- rowSums(current[,c(3:20)])
high.nn$High.nn <- rowSums(high.nn[,c(3:20)])
high.nn.q97.5$High.nn.q97.5 <- rowSums(high.nn.q97.5[,c(3:20)])
high.nn.q2.5$High.nn.q2.5 <- rowSums(high.nn.q2.5[,c(3:20)])
low.nn$Low.nn <- rowSums(low.nn[,c(3:20)])
low.nn.q97.5$Low.nn.q97.5 <- rowSums(low.nn.q97.5[,c(3:20)])
low.nn.q2.5$Low.nn.q2.5 <- rowSums(low.nn.q2.5[,c(3:20)])

birds <- data.frame(current$Atoll_Island, current$Current,
                    high.nn$High.nn, high.nn.q2.5$High.nn.q2.5, high.nn.q97.5$High.nn.q97.5,
                    low.nn$Low.nn, low.nn.q2.5$Low.nn.q2.5, low.nn.q97.5$Low.nn.q97.5)

# Tidy up:

rm(current, high.nn, low.nn, high.nn.q2.5, high.nn.q97.5, low.nn.q2.5, low.nn.q97.5)

names(birds) <- c("Atoll_Island", "Current",
                  "High.nn", "High.nn.low", "High.nn.high",
                  "Low.nn", "Low.nn.low", "Low.nn.high")

birds <- gather(birds, Estimate, Value, Current:Low.nn.high)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Split into Atoll & Island in excel

# write.csv(birds, "SNP_Data/Processed/Seabird_Numbers_Predicted_Chagos.csv")
birds <- read.csv("SNP_Data/Processed/Seabird_Numbers_Predicted_Chagos.csv")

birds$Estimate <- as.factor(birds$Estimate)
levels(birds$Estimate) <- c("Current", "Eradication", "Restoration")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Plot ####

# Convert from wide to long and facet_wrap

ggplot(birds, aes(x = Atoll_Island, )) +
  geom_bar(aes(y = Value, fill = Atoll), stat = "identity") +
  geom_errorbar(aes(ymin = Low, ymax = High, col = Atoll)) +
  scale_fill_manual(values = c("#648FFF", "#785EF0", "#DC267F", "#FE6100", "#FFB000")) +
  scale_colour_manual(values = c("#648FFF", "#785EF0", "#DC267F", "#FE6100", "#FFB000")) +
  facet_grid(Estimate~.) +
  ylab("Seabird abundance (number of pairs)") +
  xlab("Island") +
  scale_y_continuous(labels = scales::comma, trans = "log10") +
  theme_light() %+replace% theme(axis.text.x = element_text(angle = 90),
                                 panel.grid.major = element_blank(),
                                 panel.grid.minor = element_blank(),
                                 strip.text = element_text(colour = "gray30", vjust = .7),
                                 strip.background = element_rect(fill = "white", colour = "gray70"),
                                 legend.position = "bottom") +
  scale_x_discrete(labels = birds$Island)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Save:

ggsave("Plots/Supplementary/Seabird_Predictions_Chagos.png", width = 8, height = 5)
