rm(list = ls(all = TRUE))

# Data wrangling:
library(tidyr)
# For plotting:
library(ggplot2)
# (Including the bayesian stuff)
library(brms)
library(ggdist)
# (Including insets)
library(patchwork)
# (Including bringing it all together)
library(ggpubr)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Load models ####

load("SNP_ModelOutputs/Coral_Nitrogen_brms.scale.Rdata")
load("SNP_ModelOutputs/Fish_Nitrogen_brms.scale.Rdata")
load("SNP_ModelOutputs/Grazing_Nitrogen_brms.scale.Rdata")
load("SNP_ModelOutputs/Erosion_Nitrogen_brms.scale.Rdata")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Load data ####

dat <- read.csv("SNP_Data/Processed/Pred_CoralFunction_Chagos.csv")
dat <- dat[,-c(1,3:5, 7, 8, 10, 11, 13, 14, 16, 17, 19, 20, 22, 23,
               25, 26, 28, 29, 31, 32, 34, 35, 37, 38, 40, 41)]

# Wide to long
dat <- gather(dat, Scenario, Value, coralgrowth.current:erosion.gveg)
# Split up text column:
dat <- separate(dat, Scenario, into = c("Metric", "Scenario"))

# Change order of factors
dat$Scenario <- as.factor(dat$Scenario)
dat$Scenario <- factor(dat$Scenario, levels = c("current", "bveg", "gveg"))
dat$Metric <- as.factor(dat$Metric)
dat$Metric <- factor(dat$Metric, levels = c("coralgrowth", "fishbiomass", "grazing", "erosion"))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Comparison plots ####

coral.points <- ggplot(dat[which(dat$Metric=="coralgrowth"),], aes(x = Scenario, y = exp(Value))) +
  geom_jitter(alpha = 0.1, size = 2, width = 0.1) +
  # facet_wrap(.~Metric, scales = "free") +
  stat_summary(fun = "mean", aes(col = Scenario), size = .5) +
  stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 1), aes(col = Scenario),
               size = 4, alpha = 0.4, geom = "linerange") +
  theme_light() %+replace% theme(panel.grid.minor = element_blank(),
                                 legend.position = "none") +
  scale_colour_manual(values = c("#DC267F", "#785EF0", "#648FFF")) +
  scale_fill_manual(values = c("#DC267F", "#785EF0", "#648FFF")) +
  ylab(expression("Coral growth (cm"^2*" year"^-1*")")) +
  scale_x_discrete(labels = c("Current", "Rat eradication", "Rat eradication &\nhabitat restoration")) +
  guides(fill = guide_legend(title = element_blank()), col = "none") +
  xlab("Scenario")
  # expand_limits(y = c(50,180))

fish.points <- ggplot(dat[which(dat$Metric=="fishbiomass"),], aes(x = Scenario, y = exp(Value))) +
  geom_jitter(alpha = 0.1, size = 2, width = 0.1) +
  # facet_wrap(.~Metric, scales = "free") +
  stat_summary(fun = "mean", aes(col = Scenario), size = .5) +
  stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 1), aes(col = Scenario),
               size = 4, alpha = 0.4, geom = "linerange") +
  theme_light() %+replace% theme(panel.grid.minor = element_blank(),
                                 legend.position = "none") +
  scale_colour_manual(values = c("#DC267F", "#785EF0", "#648FFF")) +
  scale_fill_manual(values = c("#DC267F", "#785EF0", "#648FFF")) +
  ylab(expression("Reef fish biomass (kg ha"^-1*")")) +
  scale_x_discrete(labels = c("Current", "Rat eradication", "Rat eradication &\nhabitat restoration")) +
  guides(fill = guide_legend(title = element_blank()), col = "none") +
  xlab("Scenario")
  # expand_limits(y = c(400,1000))

graz.points <- ggplot(dat[which(dat$Metric=="grazing"),], aes(x = Scenario, y = exp(Value))) +
  geom_jitter(alpha = 0.1, size = 2, width = 0.1) +
  # facet_wrap(.~Metric, scales = "free") +
  stat_summary(fun = "mean", aes(col = Scenario), size = .5) +
  stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 1), aes(col = Scenario),
               size = 4, alpha = 0.4, geom = "linerange") +
  theme_light() %+replace% theme(panel.grid.minor = element_blank(),
                                 legend.position = "none") +
  scale_colour_manual(values = c("#DC267F", "#785EF0", "#648FFF")) +
  scale_fill_manual(values = c("#DC267F", "#785EF0", "#648FFF")) +
  ylab(expression("Grazing (% reef year"^-1*")")) +
  scale_x_discrete(labels = c("Current", "Rat eradication", "Rat eradication &\nhabitat restoration")) +
  guides(fill = guide_legend(title = element_blank()), col = "none") +
  xlab("Scenario")
  # expand_limits(y = c(4,20))

eros.points <- ggplot(dat[which(dat$Metric=="erosion"),], aes(x = Scenario, y = exp(Value))) +
  geom_jitter(alpha = 0.1, size = 2, width = 0.1) +
  # facet_wrap(.~Metric, scales = "free") +
  stat_summary(fun = "mean", aes(col = Scenario), size = .5) +
  stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 1), aes(col = Scenario),
               size = 4, alpha = 0.4, geom = "linerange") +
  theme_light() %+replace% theme(panel.grid.minor = element_blank(),
                                 legend.position = "none") +
  scale_colour_manual(values = c("#DC267F", "#785EF0", "#648FFF")) +
  scale_fill_manual(values = c("#DC267F", "#785EF0", "#648FFF")) +
  ylab(expression("Erosion (t ha"^-1*" year" ^-1*")")) +
  scale_x_discrete(labels = c("Current", "Rat eradication", "Rat eradication &\nhabitat restoration")) +
  guides(fill = guide_legend(title = element_blank()), col = "none") +
  xlab("Scenario")
  # expand_limits(y = c(0,130))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Posterior plots ####

coral.pp <- ggplot(as_draws_df(coralgrowth.model.run.scale)) +
  stat_halfeye(aes(x = b_scalelogNitrogen,  y = 0), point_interval=median_hdi,
               .width=c(.8,.5),  fill = "#FFB000", color = "#FFB000",
               fatten_point = 2, slab_alpha = .4) +
  geom_vline(xintercept = 0, lty = 2) +
  theme_light() %+replace% theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
                                 axis.ticks.y = element_blank(), axis.text.y = element_blank(),
                                 axis.title = element_blank(),
                                 # t, r, b, l
                                 plot.margin = margin(5, 7, 0, 40, unit = "pt")) +
  annotate("text", x = Inf, y = Inf, label = "Effect-size", size = 3, vjust=1.5, hjust=1.1)

coral.pp

fish.pp <- ggplot(as_draws_df(fishbiomass.model.run.scale)) +
  stat_halfeye(aes(x = b_scalelogNitrogen,  y = 0), point_interval=median_hdi,
               .width=c(.8,.5),  fill = "#FFB000", color = "#FFB000",
               fatten_point = 2, slab_alpha = .4) +
  geom_vline(xintercept = 0, lty = 2) +
  theme_light() %+replace% theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
                                 axis.ticks.y = element_blank(), axis.text.y = element_blank(),
                                 axis.title = element_blank(),
                                 plot.margin = margin(5, 7, 0, 40, unit = "pt")) +
  annotate("text", x = Inf, y = Inf, label = "Effect-size", size = 3, vjust=1.5, hjust=1.1)

fish.pp

graz.pp <- ggplot(as_draws_df(grazing.model.run.scale)) +
  stat_halfeye(aes(x = b_scalelogNitrogen,  y = 0), point_interval=median_hdi,
               .width=c(.8,.5),  fill = "#FFB000", color = "#FFB000",
               fatten_point = 2, slab_alpha = .4) +
  geom_vline(xintercept = 0, lty = 2) +
  theme_light() %+replace% theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
                                 axis.ticks.y = element_blank(), axis.text.y = element_blank(),
                                 axis.title = element_blank(),
                                 plot.margin = margin(5, 7, 0, 35, unit = "pt")) +
  annotate("text", x = Inf, y = Inf, label = "Effect-size", size = 3, vjust=1.5, hjust=1.1)


graz.pp

eros.pp <- ggplot(as_draws_df(erosion.model.run.scale)) +
  stat_halfeye(aes(x = b_scalelogNitrogen, y = .1), point_interval=median_hdi,
               .width=c(.8,.5),  fill = "#FFB000", color = "#FFB000",
               fatten_point = 2, slab_alpha = .4) +
  geom_vline(xintercept = 0, lty = 2, col = "grey70") +
  theme_light() %+replace% theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
                                 axis.ticks.y = element_blank(), axis.text.y = element_blank(),
                                 axis.title = element_blank(),
                                 plot.margin = margin(5, 7, 0, 40, unit = "pt")) +
  annotate("text", x = Inf, y = Inf, label = "Effect-size", size = 3, vjust=1.5, hjust=1.1)

eros.pp

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Save plots for manuscript ####
# l, b, r, t

# coral.plot <- coral.points + inset_element(coral.pp, 0.01, 0.6, 0.455, 0.99)
# fish.plot <- fish.points + inset_element(fish.pp, 0.01, 0.6, 0.455, 0.99)
# graz.plot <- graz.points + inset_element(graz.pp, 0.01, 0.6, 0.455, 0.99)
# eros.plot <- eros.points + inset_element(eros.pp, 0.01, 0.6, 0.455, 0.99)
# 
# ggarrange(coral.plot, fish.plot, graz.plot, eros.plot,
#           ncol = 2, nrow = 2)

ggarrange(coral.pp, fish.pp, coral.points, fish.points,
          graz.pp, eros.pp, graz.points, eros.points, 
          ncol = 2, nrow = 4,
          heights = c(0.2, 0.8, 0.2, 0.8))

ggsave("Plots/CoralFunctions_Nitrogen_Comparisons_ms.png", width = 8, height = 7)

# For a poster:

ggarrange(coral.pp, fish.pp, graz.pp, eros.pp,
          coral.points, fish.points, graz.points, eros.points, 
          ncol = 4, nrow = 2, heights = c(0.2, 0.8))

ggsave("Plots/CoralFunctions_Nitrogen_Comparisons_poster.png", width = 55, height = 15, units = "cm")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Sums for the text:

rm(list = ls(all = TRUE))

dat <- read.csv("SNP_Data/Processed/Pred_CoralFunction_Chagos.csv")
dat <- dat[,-1]

# Fish % Change:
((exp(mean(dat$fishbiomass.gveg)) - exp(mean(dat$fishbiomass.current)))/exp(mean(dat$fishbiomass.current)))*100

# Coral % change
((exp(mean(dat$coralgrowth.gveg)) - exp(mean(dat$coralgrowth.current)))/exp(mean(dat$coralgrowth.current)))*100

# Grazing % change
((exp(mean(dat$grazing.gveg)) - exp(mean(dat$grazing.current)))/exp(mean(dat$grazing.current)))*100

# Erosion % change
((exp(mean(dat$erosion.gveg)) - exp(mean(dat$erosion.current)))/exp(mean(dat$erosion.current)))*100

# Current mean fish biomass (kg ha-1) across archipelago:
exp(mean(dat$fishbiomass.current)) # kg ha-1
# Actual estiamte across entire archipelago:
exp(mean(dat$fishbiomass.current)) * 438000
# Compare with future predictions:
exp(mean(dat$fishbiomass.gveg)) * 438000
# Change
(exp(mean(dat$fishbiomass.gveg)) * 438000) - (exp(mean(dat$fishbiomass.current)) * 438000)
(exp(mean(dat$fishbiomass.gveg.low)) * 438000) - (exp(mean(dat$fishbiomass.current.low)) * 438000)
(exp(mean(dat$fishbiomass.gveg.high)) * 438000) - (exp(mean(dat$fishbiomass.current.high)) * 438000)

# Extrapolate erosion up to the entire archipelago too?
(exp(mean(dat$coralgrowth.gveg)) * 438000) - (exp(mean(dat$coralgrowth.current)) * 438000)
(exp(mean(dat$coralgrowth.gveg.low)) * 438000) - (exp(mean(dat$coralgrowth.current.low)) * 438000)
(exp(mean(dat$coralgrowth.gveg.high)) * 438000) - (exp(mean(dat$coralgrowth.current.high)) * 438000)
