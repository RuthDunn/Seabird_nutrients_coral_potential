rm(list = ls(all = TRUE))

# Packages:
library(brms)
library(plyr)
library(dplyr)
library(tidybayes)
library(ggplot2)
library(patchwork)
library(ggmcmc)
library(jtools)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Load data ####

fish <- read.csv("SNP_Data/Processed/FishBiomass_Nitrogen.csv")

# Convert all character columns to factors:
fish <- as.data.frame(unclass(fish),
                      stringsAsFactors = TRUE)

fish <- fish[complete.cases(fish), ]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Aggregate fish biomass data: ####

biomass.df <- aggregate(fish[,"Biomass"], by = list(fish$Atoll, fish$Island, fish$Transect), FUN = sum)
nitrogen.df <- aggregate(fish[,"kg_N_ha_yr"], by = list(fish$Atoll, fish$Island, fish$Transect), FUN = mean)

fish <- cbind(nitrogen.df[,c(1,2,4)], biomass.df[,4])

colnames(fish) <- c("Atoll", "Island", "Nitrogen", "FishBiomass")

rm(biomass.df, nitrogen.df)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Quick plots: ####

ggplot(data = fish, aes(x = Nitrogen, y = FishBiomass)) +
  geom_point() +
  geom_smooth(method = "lm")

ggplot(data = fish, aes(x = log(Nitrogen), y = FishBiomass)) +
  geom_point() +
  geom_smooth(method = "lm")

# ggplot(data = ch_2015_div_fish, aes(x = log_bio_kg_ha, y = log(sum_fish_bio_kg_ha))) +
#   geom_point() +
#   geom_smooth(method = "lm")

# Log both fish biomass and nitrogen input

fish$logNitrogen <- log(fish$Nitrogen)
fish$scalelogNitrogen <- scale(log(fish$Nitrogen), center = TRUE, scale = TRUE)
fish$logFishBiomass <- log(fish$FishBiomass)
fish$Atoll_Island <- paste(fish$Atoll, fish$Island, sep = "_")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Run model ####

# fishbiomass.model.run.scale <- brm(logFishBiomass ~ scalelogNitrogen + (1|Atoll_Island),
#                             data = fish,
#                             iter = 3000, warmup = 1000, chains = 4, cores = 4,
#                             control = list(adapt_delta = 0.99, max_treedepth = 12)) # Resolve divergent transitions

# fishbiomass.model.run <- brm(logFishBiomass ~ logNitrogen + (1|Atoll_Island),
#                                    data = fish,
#                                    iter = 3000, warmup = 1000, chains = 4, cores = 4,
#                                    control = list(adapt_delta = 0.99, max_treedepth = 12)) # Resolve divergent transitions

# Check it out:

# print(fishbiomass.model.run.scale)
# bayes_R2(fishbiomass.model.run.scale)

# Check diagnositcs:

# plot(fishbiomass.model.run.scale, ask = FALSE)
# pp_check(fishbiomass.model.run.scale)

# plot(fishbiomass.model.run, ask = FALSE)
# pp_check(fishbiomass.model.run)

# Cool

# save(file="SNP_ModelOutputs/Fish_Nitrogen_brms.scale.Rdata", list="fishbiomass.model.run.scale")
# save(file="SNP_ModelOutputs/Fish_Nitrogen_brms.Rdata", list="fishbiomass.model.run")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Make predictions ####

# Load models:
load("SNP_ModelOutputs/Fish_Nitrogen_brms.Rdata")
load("SNP_ModelOutputs/Fish_Nitrogen_brms.scale.Rdata")

# What proportion of the posterior distribution is above 0?
hypothesis(fishbiomass.model.run.scale, "scalelogNitrogen>0")

# Load predication data:
pred.data <- read.csv("SNP_Data/Processed/Seabird_NutrientInput_Chagos_Predicted.csv")

# High.nn scenario:
pred.data.highnn <- aggregate(pred.data$NitrogenInput.high.nn, by = list(pred.data$Atoll_Island), FUN = "sum")
names(pred.data.highnn)[1] <- "Atoll_Island"
pred.data.highnn$logNitrogen <- log(pred.data.highnn$x)

fish.highnn <- as.data.frame(predict(fishbiomass.model.run,
                                     newdata = pred.data.highnn[,c("Atoll_Island", "logNitrogen")],
                                     allow_new_levels = TRUE))

pred.data.highnn <- cbind(pred.data.highnn, fish.highnn)
rm(fish.highnn)
pred.data.highnn <- pred.data.highnn[,-2]

# Low.nn scenario:
pred.data.lownn <- aggregate(pred.data$NitrogenInput.low.nn, by = list(pred.data$Atoll_Island), FUN = "sum")
names(pred.data.lownn)[1] <- "Atoll_Island"
pred.data.lownn$logNitrogen <- log(pred.data.lownn$x)

fish.lownn <- as.data.frame(predict(fishbiomass.model.run,
                                    newdata = pred.data.lownn[,c("Atoll_Island", "logNitrogen")],
                                    allow_new_levels = TRUE))
pred.data.lownn <- cbind(pred.data.lownn, fish.lownn)
rm(fish.lownn)
pred.data.lownn <- pred.data.lownn[,-2]

write.csv(pred.data.highnn, "SNP_Data/Processed/Pred_highnn_FishBiomass_Nitrogen.csv")
write.csv(pred.data.lownn, "SNP_Data/Processed/Pred_lownn_FishBiomass_Nitrogen.csv")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Posterior plot ####

# Model variable names:
get_variables(fishbiomass.model.run.scale)

# Plot:
p2 <- ggplot(as_draws_df(fishbiomass.model.run.scale)) +
  geom_vline(xintercept = 0, lty = 2) +
  stat_halfeye(aes(x = b_scalelogNitrogen,  y = 5), point_interval=median_hdi,
               .width=c(.8,.5),  alpha = .7, fill = "#D1E1EC", fatten_point = 2, slab_alpha = .6) +
  stat_halfeye(aes(x = b_scalelogNitrogen,  y = 5), point_interval=median_hdi,
               .width=c(.8,.5),  color = "#01386B", slab_alpha = 0) +
  xlab("")+
  ylab("")+
  theme_light() %+replace% theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.ticks.y = element_blank(), axis.text.y = element_blank())


p2.t <- ggplot(as_draws_df(fishbiomass.model.run.scale)) +
  stat_halfeye(aes(y = b_scalelogNitrogen, x = 5), point_interval=median_hdi,
               .width=c(.8,.5),  alpha = .7, fill = "#c8d6fa", fatten_point = 2, slab_alpha = .4) +
  geom_hline(yintercept = 0, lty = 2) +
  stat_halfeye(aes(y = b_scalelogNitrogen, x = 5), point_interval=median_hdi,
               .width=c(.8,.5),  color = "#648FFF", slab_alpha = 0) +
  xlab("")+
  ylab("")+
  # scale_y_continuous(breaks = seq(-.4, 1.2, by = 0.4)) +
  theme_light() %+replace% theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                 axis.ticks.x = element_blank(), axis.text.x = element_blank()) +
  xlab("Effect size")

p2.t

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Conditional effects plot ####

# Model variable names:
get_variables(fishbiomass.model.run)

# Extract conditional effects:
ce.Fish.Nitrogen <- conditional_effects(fishbiomass.model.run, prob = .80,
                                 effects = 'logNitrogen', resp = "logFishBiomass",
                                 plot = FALSE)[[1]]

# Extract partialised residuals
pr.Fish.Nitrogen <- partialize(fishbiomass.model.run,
                              vars= "logNitrogen", resp = 'logFishBiomass', data = fish)

ggplot() +
  geom_ribbon(data = ce.Fish.Nitrogen, aes(x = logNitrogen, ymin = lower__, ymax=upper__),
              alpha = .4, fill = "#c8d6fa") +
  geom_jitter(data = pr.Fish.Nitrogen, aes(x = logNitrogen, y = logFishBiomass),
              alpha = 0.6, width = .1) +
  geom_point(data = pred.data.highnn, aes(x = logNitrogen, y = Estimate))

  geom_line(data = ce.Fish.Nitrogen, aes(x = logNitrogen, y = estimate__),
            color = "#648FFF") +
  xlab(expression("log Seabird nitrogen input (kg ha"^-1*" year"^-1*")")) +
  ylab(expression("log Reef fish biomass (kg ha"^-1*")")) +
  theme_bw()+
  guides(size = "none", colour = "none",  fill = guide_legend(title="Confidence level")) + 
  theme(panel.grid.major = element_blank(), # remove gridlines
        panel.grid.minor = element_blank(), #remove gridlines
        strip.background = element_blank(),
        legend.position='none')

p1

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Save plots ####

# A) Left and right?
p1 + p2.t + plot_layout(widths = c(3,1))

ggsave("Plots/Fish_Nitrogen_ppt.png", width = 5, height = 4)
