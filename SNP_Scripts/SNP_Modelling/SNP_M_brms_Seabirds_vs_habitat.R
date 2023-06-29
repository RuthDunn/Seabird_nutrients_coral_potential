
# Get ready ####

# Clear environment:

rm(list = ls(all = TRUE))

# Load packages:

library(brms)
library(dplyr)
library(tidybayes)
library(ggplot2)
library(patchwork)
library(ggmcmc)
library(jtools)
library(cmdstanr)
library(broom.mixed)
library(bayesplot)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Load & edit data ####

seabirds <- read.csv("SNP_Data/Processed/SeabirdDensityDrivers.csv")
seabirds <- seabirds[,-1]

# Remove "eradicated" and "unknown"
seabirds$Rattus_rattus <- as.factor(seabirds$Rattus_rattus)
table(seabirds$Rattus_rattus)
levels(seabirds$Rattus_rattus) <- c("A", NA, "P", NA)

# Remove rows without veg or rat data
seabirds <- seabirds[complete.cases(seabirds), ]

# Go from wide to long data format
seabirds <- gather(seabirds, Species, Pairs, Sula_sula:Anous_tenuirostris, factor_key = TRUE)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Quick plots: ####

# Island size

ggplot(data = seabirds[which(seabirds$Pairs > 0.5),], aes(x = log(Size_Ha), y = log(Pairs), fill = Rattus_rattus)) +
  geom_point() +
  geom_smooth(method = "lm")

# Veg cover

ggplot(data = seabirds[which(seabirds$Pairs > 0.5),], aes(x = (100-NonNativeForest_p), y = log(Pairs), fill = Rattus_rattus)) +
  geom_point() +
  geom_smooth(method = "lm") +
  xlab("Native forest (%)")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Transform data: ####

# Transform covariates:

seabirds$logArea <- log(seabirds$Size_Ha)
seabirds$sclogArea <- scale(log(seabirds$Size_Ha), center = TRUE, scale = TRUE)
# attr(,"scaled:center")
# [1] 2.710195 # (mean)
# attr(,"scaled:scale")
# [1] 1.832027 # (std dev)

seabirds$NativeVeg <- 100-seabirds$NonNativeForest_p
seabirds$scNativeVeg <- scale(100-seabirds$NonNativeForest_p, center = TRUE, scale = TRUE)
# attr(,"scaled:center")
# [1] 43.32374
# attr(,"scaled:scale")
# [1] 28.6727

# Data are zero-inflated:

hist(seabirds$Pairs)

hist(log(seabirds[which(seabirds$Pairs > 1),"Pairs"]))

# Use a hurdle log normal model

# https://www.andrewheiss.com/blog/2022/05/09/hurdle-lognormal-gaussian-brms/#3-hurdle-lognormal-model

# Family: hurdle_lognormal 
# Links: mu = identity; sigma = identity; hu = identity 

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Run model ####

seabirds.model.sc <- brm(bf(Pairs ~ Rattus_rattus + sclogArea + scNativeVeg + Species +
                              (1|Atoll_Island),
                            hu ~ Rattus_rattus),
                         data = seabirds,
                         family = hurdle_lognormal(),
                         iter = 3000, warmup = 1000, chains = 4, seed = 1234, silent = 2,
                         prior = c(prior(normal(6,1), class = Intercept),                           # RFB population intercept prior
                                   prior(normal(1,0.5), coef = "SpeciesOnychoprion_fuscatus"),      # Expected difference in sooty tern numbers (in comparison to RFBs)
                                   prior(normal(0,0.5), coef = "SpeciesAnous_tenuirostris"),        # Expected difference in lesser noddies numbers (in comparison to RFBs)
                                   prior(normal(-1,0.5), coef = "Rattus_rattusP"),                  # -ve influence of rats
                                   prior(normal(1,0.5), coef = "sclogArea"),                        # +ve influence of area
                                   prior(normal(1,0.5), coef = "scNativeVeg"),                      # +ve influence of native veg
                                   prior(normal(-0.4, 0.05), class = Intercept, dpar = "hu"),       # 46% of the data is 0s, so let's go with ~ 40% caused by rats
                                   prior(normal(1,0.5), coef = "Rattus_rattusP", dpar = "hu")))     # +ve influence of rats on the data being 0

plot(seabirds.model.sc)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Proportion of 0s ####

hu_intercept <- tidy(seabirds.model.sc) |> 
  filter(term == "hu_(Intercept)") |> 
  pull(estimate)

# Logit scale intercept
hu_intercept
# b_hu_Intercept 
# -1.127355 

# Transformed to a probability/proportion
plogis(hu_intercept)
# b_hu_Intercept 
# 0.2446496 

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# LOO validation
# https://mc-stan.org/loo/articles/online-only/faq.html

loo1 <- loo(seabirds.model.sc, save_psis = TRUE)

print(loo1)
plot(loo1, label_points = T)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Check model out ####

# PP Check

# https://mc-stan.org/bayesplot/articles/graphical-ppcs.html

# Exponential data:
# Logged data:
color_scheme_set("gray")
pred <- posterior_predict(seabirds.model.sc)
bayesplot::ppc_dens_overlay(y = log1p(seabirds$Pairs), 
                            yrep = log1p(pred[1:10,])) + theme_light() + ylim(c(0,0.18))

tidy(seabirds.model.sc, conf.method = "HPDinterval", conf.level = 0.95)
# Here, h_(Intercept) is the intercept for the logistic regression model
# used for the hurdle process

# Quickly check out rat effect:
# Predicted pairs (as influenced by rat status):
conditional_effects(seabirds.model.sc)
conditional_effects(seabirds.model.sc, dpar = "hu")

# What proportion of the posterior distribution is above 0?

hypothesis(seabirds.model.sc, "sclogArea>0")
hypothesis(seabirds.model.sc, "scNativeVeg>0")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Save model with variables not scaled & centred, to help with predictions ####

seabirds.model <- brm(bf(Pairs ~ Rattus_rattus + logArea + NativeVeg + Species +
                              (1|Atoll_Island),
                            hu ~ Rattus_rattus),
                         data = seabirds,
                         family = hurdle_lognormal(),
                         iter = 3000, warmup = 1000, chains = 4, seed = 1234, silent = 2,
                         prior = c(prior(normal(6,1), class = Intercept),                           # RFB population intercept prior
                                   prior(normal(1,0.5), coef = "SpeciesOnychoprion_fuscatus"),      # Expected difference in sooty tern numbers (in comparison to RFBs)
                                   prior(normal(0,0.5), coef = "SpeciesAnous_tenuirostris"),        # Expected difference in lesser noddies numbers (in comparison to RFBs)
                                   prior(normal(-1,0.5), coef = "Rattus_rattusP"),                  # -ve influence of rats
                                   prior(normal(1,0.5), coef = "logArea"),                          # +ve influence of area
                                   prior(normal(1,0.5), coef = "NativeVeg"),                        # +ve influence of native veg
                                   prior(normal(-0.4, 0.05), class = Intercept, dpar = "hu"),       # 46% of the data is 0s, so let's go with ~ 40% caused by rats
                                   prior(normal(1,0.5), coef = "Rattus_rattusP", dpar = "hu")))     # +ve influence of rats on the data being 0


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Cool, save models:

save(file="SNP_ModelOutputs/Seabirds_Habitat_brms_priors_sc.Rdata", list="seabirds.model.sc")
save(file="SNP_ModelOutputs/Seabirds_Habitat_brms_priors.Rdata", list="seabirds.model")
