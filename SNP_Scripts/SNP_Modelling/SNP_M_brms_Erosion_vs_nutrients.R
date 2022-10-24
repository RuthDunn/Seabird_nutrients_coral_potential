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

# Load data ####

fish <- read.csv("SNP_Data/Processed/CoralFunction_Nitrogen.csv")
fish <- fish[,-1]

# Convert all character columns to factors:
fish <- as.data.frame(unclass(fish),
                      stringsAsFactors = TRUE)

fish <- fish[complete.cases(fish), ]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Quick plots: ####

ggplot(data = fish, aes(x = kg_N_ha_yr, y = Erosion_t_ha)) +
  geom_point() +
  geom_smooth(method = "lm")

ggplot(data = fish, aes(x = log(kg_N_ha_yr), y = Erosion_t_ha)) +
  geom_point() +
  geom_smooth(method = "lm")

ggplot(data = fish, aes(x = log(kg_N_ha_yr), y = log(Erosion_t_ha+1))) +
  geom_point() +
  geom_smooth(method = "lm")

# Log both fish biomass and nitrogen input

fish$logNitrogen <- log(fish$kg_N_ha_yr)
fish$scalelogNitrogen <- scale(log(fish$kg_N_ha_yr), center = TRUE, scale = TRUE)
fish$log1Erosion <- log(fish$Erosion_t_ha+1)

fish$Atoll_Island <- paste (fish$Atoll, fish$Island, sep = "_")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Run model ####

# erosion.model.run.scale <- brm(log1Erosion ~ scalelogNitrogen + (1|Atoll_Island),
#                                    data = fish,
#                                    iter = 3000, warmup = 1000, chains = 4, cores = 4,
#                                    control = list(adapt_delta = 0.99, max_treedepth = 12)) # Resolve divergent transitions 

# erosion.model.run <- brm(log1Erosion ~ logNitrogen + (1|Atoll_Island),
#                              data = fish,
#                              iter = 3000, warmup = 1000, chains = 4, cores = 4,
#                              control = list(adapt_delta = 0.99, max_treedepth = 12)) # Resolve divergent transitions 

# Check it out:

# print(erosion.model.run.scale)
# print(erosion.model.run)

# Check diagnositcs:

# plot(erosion.model.run.scale, ask = FALSE)
# pp_check(erosion.model.run.scale)

# plot(erosion.model.run, ask = FALSE)
# pp_check(erosion.model.run)

# Cool

# save(file="SNP_ModelOutputs/Erosion_Nitrogen_brms.scale.Rdata", list="erosion.model.run.scale")
# save(file="SNP_ModelOutputs/Erosion_Nitrogen_brms.Rdata", list="erosion.model.run")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Make predictions ####

# Load models:
load("SNP_ModelOutputs/Erosion_Nitrogen_brms.Rdata")
load("SNP_ModelOutputs/Erosion_Nitrogen_brms.scale.Rdata")

# Load predication data:
pred.data <- read.csv("SNP_Data/Processed/Seabird_NutrientInput_Predicted.csv")

# High.nn scenario:
pred.data.highnn <- aggregate(pred.data$NitrogenInput.high.nn, by = list(pred.data$Atoll_Island), FUN = "sum")
names(pred.data.highnn)[1] <- "Atoll_Island"
pred.data.highnn$logNitrogen <- log(pred.data.highnn$x)

fish.highnn <- as.data.frame(predict(erosion.model.run,
                                     newdata = pred.data.highnn[,c("Atoll_Island", "logNitrogen")],
                                     allow_new_levels = TRUE))

pred.data.highnn <- cbind(pred.data.highnn, fish.highnn)
rm(fish.highnn)
pred.data.highnn <- pred.data.highnn[,-2]

# Low.nn scenario:
pred.data.lownn <- aggregate(pred.data$NitrogenInput.low.nn, by = list(pred.data$Atoll_Island), FUN = "sum")
names(pred.data.lownn)[1] <- "Atoll_Island"
pred.data.lownn$logNitrogen <- log(pred.data.lownn$x)

fish.lownn <- as.data.frame(predict(erosion.model.run,
                                    newdata = pred.data.lownn[,c("Atoll_Island", "logNitrogen")],
                                    allow_new_levels = TRUE))
pred.data.lownn <- cbind(pred.data.lownn, fish.lownn)
rm(fish.lownn)
pred.data.lownn <- pred.data.lownn[,-2]

write.csv(pred.data.highnn, "SNP_Data/Processed/Pred_highnn_Erosion_Nitrogen.csv")
write.csv(pred.data.lownn, "SNP_Data/Processed/Pred_lownn_Erosion_Nitrogen.csv")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Posterior plot ####

# Model variable names:
get_variables(erosion.model.run.scale)

# Plot:
p2 <- ggplot(as_draws_df(erosion.model.run.scale)) +
  stat_halfeye(aes(y = b_scalelogNitrogen,  x = 5), point_interval=median_hdi,
               .width=c(.8,.5),  alpha = .7, fill = "#f2d9a2", fatten_point = 2, slab_alpha = .4) +
  stat_halfeye(aes(y = b_scalelogNitrogen,  x = 5), point_interval=median_hdi,
               .width=c(.8,.5),  color = "#FFB000", slab_alpha = 0) +
  geom_hline(yintercept = 0, lty = 2) +
  xlab("")+
  ylab("")+
  # scale_y_continuous(breaks = seq(-1, 2, by = 0.5)) +
  theme_light() %+replace% theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                 axis.ticks.x = element_blank(), axis.text.x = element_blank())

p2

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# Conditional effects plot ####

# Model variable names:
get_variables(erosion.model.run)

# Extract conditional effects:
ce.Erosion.Nitrogen <- conditional_effects(erosion.model.run, prob = .80,
                                        effects = 'logNitrogen', resp = "sqrtErosion",
                                        plot = FALSE)[[1]]

# Extract partialised residuals
pr.Erosion.Nitrogen <- partialize(erosion.model.run,
                               vars= "logNitrogen", resp = 'sqrtErosion', data = fish)



p1 <- pr.Erosion.Nitrogen %>%
  select(logNitrogen, log1Erosion)%>%
  ggplot(aes(x = logNitrogen, y = log1Erosion)) +
  geom_ribbon(data = ce.Erosion.Nitrogen, aes(ymin = lower__, ymax=upper__),
              alpha = .4, fill = "#f2d9a2") +
  geom_jitter(alpha = 0.6, width = .1) +
  geom_line(data = ce.Erosion.Nitrogen, aes(x = logNitrogen, y = estimate__),
            color = "#FFB000") +
  xlab(expression("log Seabird nitrogen input (kg ha"^-1*" year"^-1*")")) +
  ylab(expression("log(x+1) Erosion (t ha"^-1*" year" ^-1*")")) +
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
p1 + p2 + plot_layout(widths = c(3,1))

ggsave("Plots/Erosion_Nitrogen_ppt.png", width = 5, height = 4)
