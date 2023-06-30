rm(list = ls(all = TRUE))

# Packages:
library(tidyverse)
library(brms)
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

ggplot(data = fish, aes(x = kg_N_ha_yr, y = grazing_prop_reef)) +
  geom_point() +
  geom_smooth(method = "lm")

ggplot(data = fish, aes(x = log(kg_N_ha_yr), y = grazing_prop_reef)) +
  geom_point() +
  geom_smooth(method = "lm")

ggplot(data = fish, aes(x = log(kg_N_ha_yr), y = log(grazing_prop_reef + 1))) +
  geom_point() +
  geom_smooth(method = "lm")

# Log both fish biomass and nitrogen input

fish$logNitrogen <- log(fish$kg_N_ha_yr)
fish$clogNitrogen <- scale(log(fish$kg_N_ha_yr), center = TRUE, scale = FALSE)
fish$log1Grazing <- log(fish$grazing_prop_reef+1)

fish$Atoll_Island <- paste(fish$Atoll, fish$Island, sep = "_")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Run model ####

grazing.model.run.center <- brm(log1Grazing ~ clogNitrogen + (1|Atoll_Island),
                               data = fish,
                               iter = 3000, warmup = 1000, chains = 4, cores = 4,
                               control = list(adapt_delta = 0.99, max_treedepth = 12), # Resolve divergent transitions
                               prior = c(prior(normal(2.5, 0.5), class = "Intercept"),                # Based on Graham et al, in review
                                         prior(normal(1,0.5), class = "b", coef = "clogNitrogen")))   # We expect a positive effect, like that seen with seabirds in Graham et al 2018

# Check diagnostics:

plot(grazing.model.run.center, ask = FALSE)
pp_check(grazing.model.run.center)

loo1 <- loo(grazing.model.run.center, save_psis = TRUE)

loo1
plot(loo1)

# Check model out:

tidy(grazing.model.run.center, conf.method = "HPDinterval", conf.level = 0.95)

save(file="SNP_ModelOutputs/Grazing_Nitrogen_c.Rdata", list="grazing.model.run.center")
