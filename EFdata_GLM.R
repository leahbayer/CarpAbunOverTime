# DAW/SAW Analysis Subgroup
# Upper Mississippi River Sub-basin
# EF data

# Clear all objects from the environment
rm(list=ls())

# Load required packages
library(ggplot2) # data visualization
library(MASS) # for fitting NB GLM
library(tidyverse) # data wrangling
library(emmeans) # estimating marginal means
library(pscl)


########################### Data Organization #################################
# Pull in clean data from script
source('EFdata_Clean.R')

# Convert necessary columns to factors or numeric
data$CountLarge <- as.numeric(data$CountLarge) # carp count

# Remove samples with negative or zero sampling time
data <- data %>% filter(FishTime > 0)


########################### Testing model fits ################################
# Poisson vs Negative Binomial vs Ordinary Least Squares
pois <- glm(CountLarge ~ 1 + offset(log(FishTime)), 
            family = poisson(link = "log"), data = data)
summary(pois)

nb <- glm.nb(CountLarge ~ 1 + offset(log(FishTime)), data = data)
summary(nb)


zip_mod <- zeroinfl(CountLarge ~ 1 + offset(log(FishTime)) | 1,
                    data = data,
                    dist = "poisson")

zinb_mod <- zeroinfl(CountLarge ~ 1 + offset(log(FishTime)) | 1,
                     data = data,
                     dist = "negbin")

AIC(pois, nb, zip_mod, zinb_mod)





ols <- lm(CountLarge ~ 1 + offset(log(FishTime)), data = data)
summary(ols)

AIC(pois, nb, ols) # 37715 (P) vs 15454 (NB) vs 29524 (LS)

# Dispersion
dispP <- sum(residuals(pois, type = "pearson")^2) / df.residual(pois)
dispP
# 14.3 = strong overdispersion

dispNB <- sum(residuals(nb, type = "pearson")^2) / df.residual(nb)
dispNB
# 1.1

# Plot resid vs fitted
par(mfrow = c(1, 3))
plot(fitted(pois), resid(pois, type = "pearson"), main = "Poisson Residuals")
abline(h = 0, col = "red")

plot(fitted(nb), resid(nb, type = "pearson"), main = "NB Residuals")
abline(h = 0, col = "red")

plot(fitted(ols), resid(ols), main = "OLS Residuals")
abline(h = 0, col = "red")
