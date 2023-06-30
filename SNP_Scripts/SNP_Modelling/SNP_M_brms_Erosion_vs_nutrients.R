rm(list = ls(all = TRUE))

# Packages:
library(brms)
library(dplyr)
library(tidybayes)
library(ggplot2)
library(patchwork)
library(ggmcmc)
library(jtools)
library(broom.mixed)

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

ggplot(data = fish, aes(x = log(kg_N_ha_yr), y = log(Erosion_t_ha+1))) +
  geom_point() +
  geom_smooth(method = "lm")

hist(log(fish$kg_N_ha_yr))
hist(log(fish$Erosion_t_ha+1))

# Log both fish biomass and nitrogen input
# Also center nitrogen
fish$clogNitrogen <- scale(log(fish$kg_N_ha_yr), center = TRUE, scale = FALSE)
fish$log1Erosion <- log(fish$Erosion_t_ha+1)

fish$Atoll_Island <- paste (fish$Atoll, fish$Island, sep = "_")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Run model ####

erosion.model.run.center <- brm(log1Erosion ~ clogNitrogen + (1|Atoll_Island),
                                data = fish,
                                iter = 3000, warmup = 1000, chains = 4, cores = 4,
                                control = list(adapt_delta = 0.99, max_treedepth = 12),             # Resolve divergent transitions
                                prior = c(prior(normal(4.6, 0.5), class = "Intercept"),                # Based on Graham et al 2018
                                          prior(normal(1,0.5), class = "b", coef = "clogNitrogen")))   # We expect a positive effect, like that seen with seabirds in Graham et al 2018

# Check diagnostics:

plot(erosion.model.run.center, ask = FALSE)
pp_check(erosion.model.run.center)

loo1 <- loo(erosion.model.run.center, save_psis = TRUE)

loo1
plot(loo1)

# Cool

# Check model out:

tidy(erosion.model.run.center, conf.method = "HPDinterval", conf.level = 0.95)

# What proportion of the posterior distribution is above 0?
hypothesis(erosion.model.run.center, "clogNitrogen>0")

save(file="SNP_ModelOutputs/Erosion_Nitrogen_c.Rdata", list="erosion.model.run.center")
# save(file="SNP_ModelOutputs/Erosion_Nitrogen_brms.Rdata", list="erosion.model.run")
