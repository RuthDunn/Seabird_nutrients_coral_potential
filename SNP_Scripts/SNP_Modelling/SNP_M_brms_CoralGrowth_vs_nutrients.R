rm(list = ls(all = TRUE))

# Packages:
library(tidyverse)
library(brms)
library(broom.mixed)

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
coral$clogNitrogen <- scale(log(coral$kg_N_ha_yr), center = TRUE, scale = FALSE)
# +ve values only:
coral$logCoralGrowth <- log(coral$diff_SA_year)
# (Errors are cus we're not including negative values)

coral$Atoll_Island <- paste(coral$Atoll, coral$Island, sep = "_")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Run models ####

coralgrowth.model.run.center <- brm(logCoralGrowth ~ clogNitrogen + (1|Atoll_Island),
                                    data = coral,
                                    iter = 3000, warmup = 1000, chains = 4, cores = 4,
                                    control = list(adapt_delta = 0.99, max_treedepth = 12), # Resolve divergent transitions
                                    prior = c(prior(normal(4.4, 0.5), class = "Intercept"),                # Based on Benkwitt et al, in review
                                              prior(normal(1,0.5), class = "b", coef = "clogNitrogen")))   # We expect a positive effect, like that seen with seabirds in Benkwitt et al, in review

# Check diagnostics:

plot(coralgrowth.model.run.center, ask = FALSE)
pp_check(coralgrowth.model.run.center)

loo1 <- loo(coralgrowth.model.run.center, save_psis = TRUE)

loo1
plot(loo1)

# Cool

# Check model out:

tidy(coralgrowth.model.run.center, conf.method = "HPDinterval", conf.level = 0.95)

save(file="SNP_ModelOutputs/CoralGrowth_Nitrogen_c.Rdata", list="coralgrowth.model.run.center")
