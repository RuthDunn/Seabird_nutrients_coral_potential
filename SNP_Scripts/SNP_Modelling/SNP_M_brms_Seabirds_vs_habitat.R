rm(list = ls(all = TRUE))

# Packages:
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

# Load data ####

seabirds <- read.csv("SNP_Data/Processed/SeabirdDensityDrivers.csv")
seabirds <- seabirds[,-1]

# Keep the 3 96% species:
seabirds = subset(seabirds, select = c(Atoll_Island, Rattus_rattus, Size_Ha, NonNativeForest_p,
                                       Sula_sula, Onychoprion_fuscatus, Anous_tenuirostris))

# Remove "eradicated" and "unknown"
seabirds$Rattus_rattus <- as.factor(seabirds$Rattus_rattus)
levels(seabirds$Rattus_rattus)
table(seabirds$Rattus_rattus)
levels(seabirds$Rattus_rattus) <- c("A", NA, "P", NA)

# Remove rows without veg or rat data
seabirds <- seabirds[complete.cases(seabirds), ]

# Go from wide to long
seabirds <- gather(seabirds, Species, Pairs, Sula_sula:Anous_tenuirostris, factor_key = TRUE)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Quick plots: ####

# seabirds$Pairs[seabirds$Pairs==0] <- NA
# seabirds<-seabirds[complete.cases(seabirds),]

# Island size

# ggplot(data = seabirds, aes(x = Size_Ha, y = Pairs)) +
#   geom_point() +
#   geom_smooth(method = "lm")
# 
# ggplot(data = seabirds, aes(x = log(Size_Ha), y = Pairs)) +
#   geom_point() +
#   geom_smooth(method = "lm")
# 
# ggplot(data = seabirds, aes(x = log(Size_Ha), y = log(Pairs))) +
#   geom_point() +
#   geom_smooth(method = "lm")
# 
# ggplot(data = seabirds, aes(x = log(Size_Ha), y = sqrt(Pairs))) +
#   geom_point() +
#   geom_smooth(method = "lm")
# 
# ggplot(data = seabirds, aes(x = log(Size_Ha), y = log(Pairs+1), fill = Rattus_rattus)) +
#   geom_point() +
#   geom_smooth(method = "lm")

# Veg cover

# ggplot(data = seabirds, aes(x = NonNativeForest_p, y = Pairs)) +
#   geom_point() +
#   geom_smooth(method = "lm")
# 
# ggplot(data = seabirds, aes(x = log(NonNativeForest_p), y = Pairs)) +
#   geom_point() +
#   geom_smooth(method = "lm")
# 
# ggplot(data = seabirds, aes(x = log(NonNativeForest_p), y = log(Pairs))) +
#   geom_point() +
#   geom_smooth(method = "lm")
# 
# ggplot(data = seabirds, aes(x = log(NonNativeForest_p), y = sqrt(Pairs))) +
#   geom_point() +
#   geom_smooth(method = "lm")
# 
# ggplot(data = seabirds, aes(x = log(NonNativeForest_p), y = log(Pairs+1))) +
#   geom_point() +
#   geom_smooth(method = "lm")
# 
# ggplot(data = seabirds, aes(x = log(100-NonNativeForest_p), y = log(Pairs+1), fill = Rattus_rattus)) +
#   geom_point() +
#   geom_smooth(method = "lm")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Transform data: ####

# Transform covariates:

seabirds$logArea <- log(seabirds$Size_Ha)
seabirds$sclogArea <- scale(log(seabirds$Size_Ha), center = TRUE, scale = TRUE)

seabirds$nonnonVeg <- 100-seabirds$NonNativeForest_p
seabirds$logVeg <- log(100-seabirds$NonNativeForest_p)
seabirds$sclogVeg <- scale(log(100-seabirds$NonNativeForest_p), center = TRUE, scale = TRUE)

# Transform response:

# Check for zero-inflated data:
hist(seabirds$Pairs)

# Check for over-dispersion:
mean(seabirds$Pairs)
var(seabirds$Pairs) 
# There is substantial over-dispersion

# Do other transformations:
seabirds$logPairs <- log(seabirds$Pairs+1)
seabirds$Pairs.p1 <- seabirds$Pairs+1

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Run model ####

seabirds.model.run.scale11 <- brm(bf(Pairs ~ Rattus_rattus + sclogArea + sclogVeg + Species +
                                       (1|Atoll_Island),
                                     hu ~ Rattus_rattus),
                                 data = seabirds,
                                 family = hurdle_lognormal(),
                                 iter = 3000, warmup = 1000, chains = 4, seed = 1234, silent = 2,
                                 save_all_pars = T)

# Model comparison:
# loo9 <- loo(seabirds.model.run.scale9, save_psis = TRUE)
# loo_compare(loo5, loo6, loo7, loo8, loo9)
# (Top one = best one)
# 9 is the best one:
# Family: hurdle_lognormal 
# Links: mu = identity; sigma = identity; hu = identity 

# https://www.andrewheiss.com/blog/2022/05/09/hurdle-lognormal-gaussian-brms/#3-hurdle-lognormal-model

# Run with cmdstanr (faster + more modern than rstan)

# seabirds.model.run11 <- brm(bf(Pairs ~ Rattus_rattus + logArea + logVeg + Species +
#                                  (1|Atoll_Island),
#                                hu ~ Rattus_rattus),
#                             data = seabirds,
#                             family = hurdle_lognormal(),
#                             chains = 4, iter = 3000, warmup = 1000, seed = 1234,  silent = 2)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

load("SNP_ModelOutputs/Seabirds_Habitat_brms.hln.Rdata")
load("SNP_ModelOutputs/Seabirds_Habitat_brms.scale.hln.Rdata")

tidy(seabirds.model.run11)
# Here, h_(Intercept) is the intercept for the logistic regression model
# used for the hurdle process

# Get proportion of zeros:
hu_intercept <- tidy(seabirds.model.run11) |> 
  filter(term == "hu_(Intercept)") |> 
  pull(estimate)
plogis(hu_intercept)

# Confirm this:
seabirds$is_zero <- ifelse(seabirds$Pairs == 0, TRUE, FALSE)
seabirds |> 
  count(is_zero) |> 
  mutate(prop = n / sum(n))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

conditional_effects(seabirds.model.run11, effects = "Rattus_rattus")
conditional_effects(seabirds.model.run11, dpar = "hu")

kfold1 <- kfold(seabirds.model.run11, chains = 1)
kfold1
plot(kfold1$pointwise)
# https://vasishth.github.io/bayescogsci/book/ch-cv.html

# PP Check
# Exponential data:
pp_check(seabirds.model.run11)
# Logged data:
pred <- posterior_predict(seabirds.model.run.scale1)
bayesplot::ppc_dens_overlay(y = log1p(seabirds$Pairs), 
                            yrep = log1p(pred[1:10,]))

tidy(seabirds.model.run11)

# Predicted pairs (as influenced by rat status):
conditional_effects(seabirds.model.run11, effects = "Rattus_rattus")
# Predicted probability of seeing no pairs (as influenced by rat status):
conditional_effects(seabirds.model.run11, dpar = "hu")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Cool

# save(file="SNP_ModelOutputs/Seabirds_Habitat_brms.scale.hln.Rdata", list="seabirds.model.run.scale11")
# save(file="SNP_ModelOutputs/Seabirds_Habitat_brms.hln.Rdata", list="seabirds.model.run11")

