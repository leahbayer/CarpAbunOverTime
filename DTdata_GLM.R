# DAW/SAW Analysis Subgroup
# Upper Mississippi River Sub-basin
# Ask 2

# Clear all objects from the environment
rm(list=ls())

library(MASS) # for fitting NB GLM
library(parallel) # parallel computing
library(tidyverse) # data wrangling
library(tools)

##################### Data Organization & Model Fit ###########################
# Pull in clean data from script
source('DTdata_Clean.R')

# separate analyses for high and medium abundance
# variance might be different between the two groups
data_high <- data %>% filter(Abundance == 'High')
data_mod <- data %>% filter(Abundance == 'Moderate')
data_rare <- data %>% filter(Abundance == 'Rare')

# Fit GLM
nb_high <- glm.nb(CountLarge ~ 1 + offset(log(FishTime)), data = data_high)
summary(nb_high)

nb_mod <- glm.nb(CountLarge ~ 1 + offset(log(FishTime)), data = data_mod)
summary(nb_mod)

nb_rare <- glm.nb(CountLarge ~ 1 + offset(log(FishTime)), data = data_rare)
summary(nb_rare)


########################### Testing model fits ################################
# Poisson vs Negative Binomial vs Ordinary Least Squares
pois_high <- glm(CountLarge ~ 1 + offset(log(FishTime)), 
                family = poisson(link = "log"), data = data_high)

nb_high <- glm.nb(CountLarge ~ 1 + offset(log(FishTime)), data = data_high)

AIC(pois_high, nb_high)


pois_mod <- glm(CountLarge ~ 1 + offset(log(FishTime)), 
                 family = poisson(link = "log"), data = data_mod)

nb_mod <- glm.nb(CountLarge ~ 1 + offset(log(FishTime)), data = data_mod)

AIC(pois_mod, nb_mod)


pois_rare <- glm(CountLarge ~ 1 + offset(log(FishTime)), 
                 family = poisson(link = "log"), data = data_rare)

nb_rare <- glm.nb(CountLarge ~ 1 + offset(log(FishTime)), data = data_rare)

AIC(pois_rare, nb_rare)


library(pscl)

# Zero-inflated Poisson
zip_high <- zeroinfl(CountLarge ~ 1 + offset(log(FishTime)) | 1,
                     data = data_high,
                     dist = "poisson")

# Zero-inflated Negative Binomial
zinb_high <- zeroinfl(CountLarge ~ 1 + offset(log(FishTime)) | 1,
                      data = data_high,
                      dist = "negbin")

AIC(pois_high, nb_high, zip_high, zinb_high)

zip_mod <- zeroinfl(CountLarge ~ 1 + offset(log(FishTime)) | 1,
                    data = data_mod,
                    dist = "poisson")

zinb_mod <- zeroinfl(CountLarge ~ 1 + offset(log(FishTime)) | 1,
                     data = data_mod,
                     dist = "negbin")

AIC(pois_mod, nb_mod, zip_mod, zinb_mod)

zip_rare <- zeroinfl(CountLarge ~ 1 + offset(log(FishTime)) | 1,
                     data = data_rare,
                     dist = "poisson")

zinb_rare <- zeroinfl(CountLarge ~ 1 + offset(log(FishTime)) | 1,
                      data = data_rare,
                      dist = "negbin")

AIC(pois_rare, nb_rare, zip_rare, zinb_rare)



#ols_high <- lm(CountLarge ~ Abundance + offset(log(FishTime)), data = data)

#AIC(pois, nb, ols) # 39716 (P) vs 16272 (NB) vs 22481 (LS)

# Dispersion
dispP <- sum(residuals(pois, type = "pearson")^2) / df.residual(pois)
dispP
# 17.4 = strong overdispersion

dispNB <- sum(residuals(nb, type = "pearson")^2) / df.residual(nb)
dispNB
# 1.2

# Plot resid vs fitted
par(mfrow = c(1, 3))
plot(fitted(pois), resid(pois, type = "pearson"), main = "Poisson Residuals")
abline(h = 0, col = "red")

plot(fitted(nb), resid(nb, type = "pearson"), main = "NB Residuals")
abline(h = 0, col = "red")

plot(fitted(ols), resid(ols), main = "OLS Residuals")
abline(h = 0, col = "red")





abun_emm <- emmeans(nb, ~ Abundance, offset = log(1))
abun_summary <- as.data.frame(abun_emm)

abun_summary <- abun_summary %>%
  mutate(
    y = exp(emmean) * 60,
    low = exp(asymp.LCL) * 60,
    high = exp(asymp.UCL) * 60
  )

ggplot(abun_summary, aes(x = Abundance, y = y)) +
  geom_point() +
  geom_errorbar(aes(ymin = low, ymax = high), width = 0.2) +
  labs(title = "Predicted CPUE by Abundance Category in UMR", x = NULL, 
       y = "CPUE (catch per hour)") +
  theme(
    axis.text.x = element_text(size = 20, face = "bold"),
    axis.text.y = element_text(size = 20),
    axis.title.y = element_text(size = 22, face = "bold"),
    axis.title.x = element_text(size = 22, face = "bold"),
    plot.title = element_text(size = 24, hjust = 0.5),
    strip.text = element_text(size = 20, face = "bold"),
    panel.spacing = unit(3, "lines")
  )







# new
get_preds <- function(model, abun_label) {
  emm <- emmeans(model, ~1, offset = log(1)) |> as.data.frame()
  emm |> 
    mutate(
      Abundance = abun_label,
      y    = exp(emmean) * 60,
      low  = exp(asymp.LCL) * 60,
      high = exp(asymp.UCL) * 60
    ) |> 
    dplyr::select(Abundance, y, low, high)
}

abun_summary <- bind_rows(
  get_preds(nb_rare, "Rare"),
  get_preds(nb_mod, "Median"),
  get_preds(nb_high, "High")
)

# Plot
ggplot(abun_summary, aes(x = Abundance, y = y)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = low, ymax = high), width = 0.2) +
  labs(
    title = "Predicted CPUE by Abundance Category in UMR",
    x = NULL, 
    y = "CPUE (catch per hour)"
  ) +
  theme(
    axis.text.x  = element_text(size = 20, face = "bold"),
    axis.text.y  = element_text(size = 20),
    axis.title.y = element_text(size = 22, face = "bold"),
    axis.title.x = element_text(size = 22, face = "bold"),
    plot.title   = element_text(size = 24, hjust = 0.5),
    strip.text   = element_text(size = 20, face = "bold"),
    panel.spacing = unit(3, "lines")
  )
