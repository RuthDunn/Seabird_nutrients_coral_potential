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

coral <- read.csv("SNP_Data/Processed/CoralGrowth_Nitrogen.csv")

# Convert all character columns to factors:
coral <- as.data.frame(unclass(coral),
                      stringsAsFactors = TRUE)

coral <- coral[complete.cases(coral), ]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Quick plots: ####

ggplot(data = coral, aes(x = kg_N_ha_yr, y = diff_SA_year)) +
  geom_point() +
  geom_smooth(method = "lm")

ggplot(data = coral, aes(x = log(kg_N_ha_yr), y = diff_SA_year)) +
  geom_point() +
  geom_smooth(method = "lm")

ggplot(data = coral, aes(x = log(kg_N_ha_yr), y = log(diff_SA_year + 1))) +
  geom_point() +
  geom_smooth(method = "lm")

# Log both coral growth and nitrogen input

coral$logNitrogen <- log(coral$kg_N_ha_yr)
coral$scalelogNitrogen <- scale(log(coral$kg_N_ha_yr), center = TRUE, scale = TRUE)
# +ve values only:
coral$logCoralGrowth <- log(coral$diff_SA_year)
# (Errors are cus we're not including negative values)

coral$Atoll_Island <- paste(coral$Atoll, coral$Island, sep = "_")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Run models ####

# coralgrowth.model.run.scale <- brm(logCoralGrowth ~ scalelogNitrogen + (1|Atoll_Island),
#                                    data = coral,
#                                    iter = 3000, warmup = 1000, chains = 4, cores = 4,
#                                    control = list(adapt_delta = 0.99, max_treedepth = 12)) # Resolve divergent transitions 

# coralgrowth.model.run <- brm(logCoralGrowth ~ logNitrogen + (1|Atoll_Island),
#                                    data = coral,
#                                    iter = 3000, warmup = 1000, chains = 4, cores = 4,
#                                    control = list(adapt_delta = 0.99, max_treedepth = 12)) # Resolve divergent transitions 

# Check it out:

# print(coralgrowth.model.run.scale)
# print(coralgrowth.model.run)

# Check diagnositcs:

# plot(coralgrowth.model.run.scale, ask = FALSE)
# pp_check(coralgrowth.model.run.scale)

# plot(coralgrowth.model.run, ask = FALSE)
# pp_check(coralgrowth.model.run)

# Cool

# save(file="SNP_ModelOutputs/Coral_Nitrogen_brms.scale.Rdata", list="coralgrowth.model.run.scale")
# save(file="SNP_ModelOutputs/Coral_Nitrogen_brms.Rdata", list="coralgrowth.model.run")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Make predictions ####

load("SNP_ModelOutputs/Coral_Nitrogen_brms.scale.Rdata")
load("SNP_ModelOutputs/Coral_Nitrogen_brms.Rdata")

# Load predication data:
pred.data <- read.csv("SNP_Data/Processed/Seabird_NutrientInput_Chagos_Predicted.csv")

# High.nn scenario:
pred.data.highnn <- aggregate(pred.data$NitrogenInput.high.nn, by = list(pred.data$Atoll_Island), FUN = "sum")
names(pred.data.highnn)[1] <- "Atoll_Island"
pred.data.highnn$logNitrogen <- log(pred.data.highnn$x)

coral.highnn <- as.data.frame(predict(coralgrowth.model.run,
                                     newdata = pred.data.highnn[,c("Atoll_Island", "logNitrogen")],
                                     allow_new_levels = TRUE))

pred.data.highnn <- cbind(pred.data.highnn, coral.highnn)
rm(coral.highnn)
pred.data.highnn <- pred.data.highnn[,-2]

# Low.nn scenario:
pred.data.lownn <- aggregate(pred.data$NitrogenInput.low.nn, by = list(pred.data$Atoll_Island), FUN = "sum")
names(pred.data.lownn)[1] <- "Atoll_Island"
pred.data.lownn$logNitrogen <- log(pred.data.lownn$x)

coral.lownn <- as.data.frame(predict(coralgrowth.model.run,
                                    newdata = pred.data.lownn[,c("Atoll_Island", "logNitrogen")],
                                    allow_new_levels = TRUE))
pred.data.lownn <- cbind(pred.data.lownn, coral.lownn)
rm(coral.lownn)
pred.data.lownn <- pred.data.lownn[,-2]

write.csv(pred.data.highnn, "SNP_Data/Processed/Pred_highnn_CoralGrowth_Nitrogen.csv")
write.csv(pred.data.lownn, "SNP_Data/Processed/Pred_lownn_CoralGrowth_Nitrogen.csv")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Posterior plot ####

# Model variable names:
get_variables(coralgrowth.model.run.scale)

p2 <- ggplot(as_draws_df(coralgrowth.model.run.scale)) +
  geom_vline(xintercept = 0, lty = 2) +
  stat_halfeye(aes(x = b_scalelogNitrogen,  y = 5), point_interval=median_hdi,
               .width=c(.8,.5),  alpha = .7, fill = "#D1E1EC", fatten_point = 2, slab_alpha = .6) +
  stat_halfeye(aes(x = b_scalelogNitrogen,  y = 5), point_interval=median_hdi,
               .width=c(.8,.5),  color = "#01386B", slab_alpha = 0) +
  xlab("")+
  ylab("")+
  theme_light() %+replace% theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                 axis.ticks.y = element_blank(), axis.text.y = element_blank())

p2

p2.t <- ggplot(as_draws_df(coralgrowth.model.run.scale)) +
  stat_halfeye(aes(y = b_scalelogNitrogen, x = 5), point_interval=median_hdi,
               .width=c(.8,.5),  alpha = .7, fill = "#fad0b6", fatten_point = 2, slab_alpha = .4) +
  geom_hline(yintercept = 0, lty = 2) +
  stat_halfeye(aes(y = b_scalelogNitrogen, x = 5), point_interval=median_hdi,
               .width=c(.8,.5),  color = "#FE6100", slab_alpha = 0) +
  xlab("") +
  ylab("") +
  scale_y_continuous(limits = c(-1, 1.4), breaks = seq(-1, 1.4, by = 0.4)) +
  theme_light() %+replace% theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                 axis.ticks.x = element_blank(), axis.text.x = element_blank())

p2.t

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Conditional effects plot

# Extract conditional effects:
ce.Coral.Nitrogen <- conditional_effects(coralgrowth.model.run, prob = .80,
                                        effects = 'logNitrogen', resp = "logCoralGrowth",
                                        plot = FALSE)[[1]]

# Extract partialised residuals
pr.Coral.Nitrogen <- partialize(coralgrowth.model.run,
                               vars= "logNitrogen", resp = 'logCoralGrowth', data = coral)



p1 <- pr.Coral.Nitrogen %>%
  select(logNitrogen, logCoralGrowth)%>%
  ggplot(aes(x = logNitrogen, y = logCoralGrowth)) +
  geom_ribbon(data = ce.Coral.Nitrogen, aes(ymin = lower__, ymax=upper__),
              alpha = .4, fill = "#fad0b6") +
  geom_jitter(alpha = 0.6, width = .1) +
  geom_line(data = ce.Coral.Nitrogen, aes(x = logNitrogen, y = estimate__),
            color = "#FE6100") +
  xlab(expression("log Seabird nitrogen input (kg ha"^-1*" year"^-1*")")) +
  ylab(expression("log Annual coral growth (cm"^2*" year"^-1*")")) +
  guides(size = "none", colour = "none",  fill = guide_legend(title="Confidence level")) +
  theme_light() %+replace% theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

p1
# ggsave("Plots/Coral_Nitrogen_ppt2.png", width = 5, height = 4)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Save plots ####

# A) Left and right?
p1 + p2.t + plot_layout(widths = c(3,1))

ggsave("Plots/Coral_Nitrogen_ppt.png", width = 5, height = 4)
