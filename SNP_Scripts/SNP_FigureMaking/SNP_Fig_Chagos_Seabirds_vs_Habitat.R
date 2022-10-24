rm(list = ls(all = TRUE))

# Packages:
library(brms)
library(ggplot2)
library(patchwork)
library(tidyr)
library(ggdist)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Load data ####

seabirds <- read.csv("SNP_Data/Processed/SeabirdDensityDrivers.csv")
seabirds <- seabirds[,-1]

# Keep the 3 96% species:
seabirds = subset(seabirds, select = c(Atoll_Island, Rattus_rattus, Size_Ha, NonNativeForest_p,
                                       Sula_sula, Onychoprion_fuscatus, Anous_tenuirostris))

# Change "eradicated" to "absent"
seabirds$Rattus_rattus <- as.factor(seabirds$Rattus_rattus)
levels(seabirds$Rattus_rattus) <- c("A", "A", "P", NA)

# Remove rows without veg or rat data
seabirds <- seabirds[complete.cases(seabirds), ]

# Go from wide to long
seabirds <- gather(seabirds, Species, Pairs, Sula_sula:Anous_tenuirostris, factor_key = TRUE)

# Transform data
seabirds$logPairs <- log(seabirds$Pairs+1)
seabirds$logArea <- log(seabirds$Size_Ha)
seabirds$logVeg <- log(100-seabirds$NonNativeForest_p)

# Load models
load("SNP_ModelOutputs/Seabirds_Habitat_brms.hln.Rdata")
load("SNP_ModelOutputs/Seabirds_Habitat_brms.scale.hln.Rdata")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Posterior plot ####

# Model variable names:
# get_variables(seabirds.model.run.scale)

# Plot:
p2 <- ggplot(as_draws_df(seabirds.model.run.scale11)) +
  geom_hline(yintercept = 0, lty = 2) +
  
  # Posteriors
  stat_halfeye(aes(y = b_Rattus_rattusP,  x = -1), point_interval=median_hdi,
               .width=c(.8,.5),  alpha = .7, fill = "#DC267F", fatten_point = 2, slab_alpha = .3) +
  stat_halfeye(aes(y = b_sclogArea,  x = -2), point_interval=median_hdi,
               .width=c(.8,.5),  alpha = .7, fill = "#785EF0", fatten_point = 2, slab_alpha = .2) +
  stat_halfeye(aes(y = b_sclogVeg,  x = -3), point_interval=median_hdi,
               .width=c(.8,.5),  alpha = .7, fill = "#648FFF", fatten_point = 2, slab_alpha = .2) +
  # Points
  stat_halfeye(aes(y = b_Rattus_rattusP,  x = -1), point_interval=median_hdi,
               .width=c(.8,.5),  color = "#DC267F", slab_alpha = 0) +
  stat_halfeye(aes(y = b_sclogArea,  x = -2), point_interval=median_hdi,
               .width=c(.8,.5),  color = "#785EF0", slab_alpha = 0) +
  stat_halfeye(aes(y = b_sclogVeg,  x = -3), point_interval=median_hdi,
               .width=c(.8,.5),  color = "#648FFF", slab_alpha = 0) +
  
  xlab("Scenario")+
  ylab("Effect-size")+
  theme_light() %+replace% theme(panel.grid.minor = element_blank(),
                                 axis.ticks.y = element_blank(),
                                 axis.text.x = element_text(angle = 270)) +
  scale_x_continuous(breaks = c(-1,-2,-3), labels = c("", "", ""))

p2

# Posterior plot only:

p3 <- ggplot(as_draws_df(seabirds.model.run.scale11)) +
  geom_vline(xintercept = 0, lty = 2) +
  
  # Posteriors
  stat_halfeye(aes(x = b_Rattus_rattusP,  y = -1), point_interval=median_hdi,
               .width=c(.8,.5),  alpha = .7, fill = "#DC267F", fatten_point = 2, slab_alpha = .3) +
  stat_halfeye(aes(x = b_sclogArea,  y = -2), point_interval=median_hdi,
               .width=c(.8,.5),  alpha = .7, fill = "#785EF0", fatten_point = 2, slab_alpha = .2) +
  stat_halfeye(aes(x = b_sclogVeg,  y = -3), point_interval=median_hdi,
               .width=c(.8,.5),  alpha = .7, fill = "#648FFF", fatten_point = 2, slab_alpha = .2) +
  # Points
  stat_halfeye(aes(x = b_Rattus_rattusP,  y = -1), point_interval=median_hdi,
               .width=c(.8,.5),  color = "#DC267F", slab_alpha = 0) +
  stat_halfeye(aes(x = b_sclogArea,  y = -2), point_interval=median_hdi,
               .width=c(.8,.5),  color = "#785EF0", slab_alpha = 0) +
  stat_halfeye(aes(x = b_sclogVeg,  y = -3), point_interval=median_hdi,
               .width=c(.8,.5),  color = "#648FFF", slab_alpha = 0) +
  
  ylab("Influences of seabird abundance") +
  xlab("Effect-size") +
  theme_light() %+replace% theme(panel.grid.minor = element_blank(),
                                 axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0), angle = 90)) +
  scale_y_continuous(breaks = c(-1,-2,-3), labels = c("Rat presence", "Island area", "Native habitat cover"))

p3

ggsave("Plots/Seabirds_Habitat_ms.png", width = 8, height = 3)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Conditional effects plot ####

# Model variable names:
# get_variables(seabirds.model.run)

# Plot

# plot(conditional_effects(seabirds.model.run11, effects = "Rattus_rattus"))
# plot(conditional_effects(seabirds.model.run11, effects = "logArea"))
# plot(conditional_effects(seabirds.model.run11, effects = "logVeg"))
# plot(conditional_effects(seabirds.model.run11, effects = "Species"))

p <- plot(conditional_effects(seabirds.model.run11,
                              effects = "logArea:logVeg"),
          points = T,
          point_args = c(alpha = 2/3, size = 1), mean = F)

# Data frame of conditional effects:
ce.Seabird.model <- p$`logArea:logVeg`$data

# Create a new sequence of predictor values:
ce.Seabird.model$Area_new <- seq(min(seabirds$logArea), max(seabirds$logArea), length.out = nrow(ce.Seabird.model))

p1 <- ggplot(data = seabirds, aes(x = exp(logArea))) +
  # Confidence intervals:
  geom_ribbon(data = ce.Seabird.model, aes(ymin = (lower__), ymax = (upper__), fill = effect2__), alpha = .2) +
  scale_fill_manual(values = c("#FE6100", "#FFB000", "#648FFF"), labels = c("High", "Medium", "Low")) +
  #Raw data:
  geom_point(aes(y = exp(logPairs), shape = Rattus_rattus), alpha = 0.5, size = 2) +
  scale_shape_manual(name = "", values = c(1, 16), labels = c("Absent", "Present")) +
  # Model fit lines:
  geom_line(data = ce.Seabird.model, aes(x = exp(Area_new), y = (estimate__), col = effect2__)) +
  scale_colour_manual(values = c("#FE6100", "#FFB000", "#648FFF")) +
  # Theme stuff:
  xlab("Island area (ha)") +
  ylab("Seabird abundance (number of pairs)") +
  theme_bw()+
  guides(size = "none", colour = "none",  fill = guide_legend(title="Confidence level")) + 
  theme_light() %+replace% theme(panel.grid.minor = element_blank(), #remove gridlines
        strip.background = element_blank(),
        legend.position='right') +
  guides(fill = guide_legend(order = 2, title = paste("Proportion of non-\nplantation forest")),
         shape = guide_legend(order = 1, title = "Rat status")) +
  scale_x_continuous(trans = "log10") +
  scale_y_continuous(labels = scales::comma, trans = "log10")
  
p1

# This makes sense: the slope is higher, the higher the proportion of NON non-native vegetation.

# ggsave("Plots/Seabirds_Habitat_ms.png", width = 8, height = 5)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Save plot for manuscript ####

# A) Side by side?
p2 + p1 +
  plot_layout(widths = c(1,2))

ggsave("Plots/Seabirds_Habitat_ms.png", width = 8, height = 5)
# ggsave("Plots/Seabirds_Habitat_ppt.png", width = 5, height = 4)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Values for manuscript text: ####

print(seabirds.model.run.scale11)

