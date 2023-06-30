rm(list = ls(all = TRUE))

# Packages:
library(tidyverse)
library(brms)
library(broom.mixed)

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

# Log both fish biomass and nitrogen input

ggplot(data = fish, aes(x = log(Nitrogen), y = log(FishBiomass))) +
  geom_point() +
  geom_smooth(method = "lm")

fish$logNitrogen <- log(fish$Nitrogen)
fish$clogNitrogen <- scale(log(fish$Nitrogen), center = TRUE, scale = FALSE)
fish$logFishBiomass <- log(fish$FishBiomass)
fish$Atoll_Island <- paste(fish$Atoll, fish$Island, sep = "_")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Run models ####

fishbiomass.model.run.center <- brm(logFishBiomass ~ clogNitrogen + (1|Atoll_Island),
                            data = fish,
                            iter = 3000, warmup = 1000, chains = 4, cores = 4,
                            control = list(adapt_delta = 0.99, max_treedepth = 12), # Resolve divergent transitions
                            prior = c(prior(normal(6.5, 0.5), class = "Intercept"),                # Based on Graham et al, in review
                                      prior(normal(1,0.5), class = "b", coef = "clogNitrogen")))   # We expect a positive effect, like that seen with seabirds in Graham et al 2018


# Check diagnostics:

plot(fishbiomass.model.run.center, ask = FALSE)
pp_check(fishbiomass.model.run.center)

loo1 <- loo(fishbiomass.model.run.center, save_psis = TRUE)

loo1
plot(loo1)

# Check model out:

tidy(fishbiomass.model.run.center, conf.method = "HPDinterval", conf.level = 0.95)

save(file="SNP_ModelOutputs/Fish_Nitrogen_c.Rdata", list="fishbiomass.model.run.center")
