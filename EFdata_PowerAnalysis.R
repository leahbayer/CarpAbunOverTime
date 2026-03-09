# DAW/SAW Analysis Subgroup
# Power Analysis for Upper Mississippi River Sub-basin
# Sampling effort needed for:
#     60, 70, & 80% power 
#     to detect a - 15, 25, & 50% change in 
#     adult Silver Carp abundance over 
#     a 3- and 5-year time span at 
#     alpha = 0.05
# Using EF data

# Clear all objects from the environment
rm(list=ls())

library(MASS) # for fitting NB GLM
library(parallel) # parallel computing
library(tidyverse) # data wrangling
library(tools)


##################### Data Organization & Model Fit ###########################
# Pull in clean data from script
source('EFdata_Clean.R')

# Convert necessary columns to factors or numeric
data$CountLarge <- as.numeric(data$CountLarge) # carp count

# Remove samples with negative or zero sampling time
data <- data %>% filter(FishTime > 0)

# Fit GLM
nb <- glm.nb(CountLarge ~ 1 + offset(log(FishTime)), data = data)
summary(nb)


######################## Define Functions #####################################

# Simulate count data from NB model
# applies offset & effect size to get new expected count
# returns data frame with counts over time

simulate_counts <- function(n, params, eff, theta, offset = log(15),
                            time = 1:3) {
  # n - sample size
  # params - coefficients of nb_all model
  # eff - effect size (influence of params on outcome)
  # theta - dispersion from model
  # offset - effort offset (log FishTime), multiplies by standard sampling
  # defaults to 15 minutes for electrofishing data
  # time - time effect
  
  # Expected count
  counts <- matrix(NA, nrow = n, ncol = length(time))
  
  for (yr in 1:length(time)) {
    t <- time[yr]
    exp_count <- exp(params + offset + eff * t)
    counts[, yr] <- rnbinom(n, mu = exp_count, size = theta)
  }
  
  sim_data <- data.frame(y = as.vector(counts), Year = rep(time, each = n),
                         offset = offset)
  return(sim_data)
}


# Simulate n_sim datasets for a range of sample sizes
# fits GLM to each simulated dataset
# checks whether true trend is detected (p < 0.05)
# returns power estimates for each combo

simulated_power <- function(sample_size, n_sim = 1000, params, eff, theta,
                            alpha = 0.05, offset = log(15), time){
  # sample_size - annual sample size per simulation
  # n_sim - number of simulations, defaults to 1000
  # params - intercept of log-linear model
  # eff - slope coefficient of log-linear model
  # theta - negative binomial dispersion parameter
  # alpha - type - 1 error, defaults to 0.05
  # offset - effort offset (log FishTime), multiplies by standard sampling
  # defaults to 5 minutes for dozer trawl data
  # time - vector, which years sampling occurs
  
  p_val <- numeric(n_sim) # for storing simulated p-value
  
  for(i in 1:n_sim){
    count <- simulate_counts(n = sample_size,
                             params = params,
                             eff = eff,
                             theta = theta,
                             offset = offset,
                             time = time)
    
    # testing for warnings or errors from the fitted model
    test <- try(assertCondition(glm.nb(y ~ Year + offset(offset),
                                       data = count)), silent = T)
    
    # simulating data until no warnings
    while(any(class(test[[1]]) %in%
              c('simpleWarning', 'simpleError'))){
      count <- simulate_counts(n = sample_size,
                               params = params,
                               eff = eff,
                               theta = theta,
                               offset = offset,
                               time = time)
      
      # testing for warnings or errors from the fitted model
      test <- try(assertCondition(glm.nb(y ~ Year + offset(offset),
                                         data = count)), silent = T)
    }
    
    fit <- glm.nb(y ~ Year + offset(offset), data = count)
    
    p_val[i] <- coef(summary(fit))["Year", "Pr(>|z|)"]
  }
  
  return(mean(p_val < alpha))
  
}


##################### Fit models to determine power ###########################

# Assign params
samples <- seq(25, 1000, by = 5)

# defining simulation scenarios
trends <- c(log(0.85), log(0.75), log(0.50))
names(trends) <- c('15% Reduction', '25% Reduction', '50% Reduction')

years_sampled <- list(1:3, 1:5, c(1, 3, 5))
names(years_sampled) <- c('Annually for 3 years', 'Annually for 5 years',
                          'Year 1, 3, and 5')

coef_levels <- coef(nb)
names(coef_levels) <- 'Median Abundance'

theta_levels <- nb$theta
names(theta_levels) <- 'Median Abundance'

# creating a data.frame to store output
out <- expand.grid(Trend = names(trends), Sampling = names(years_sampled),
                   Abundance = names(coef_levels)) %>%
  add_column(`Sample Size, Power = 0.60` = NA,
             `Sample Size, Power = 0.70` = NA,
             `Sample Size, Power = 0.80` = NA)

for(tr in 1:length(trends)){
  for(ys in 1:length(years_sampled)){
    for(ab in 1:length(coef_levels)){
      
      c1 <- makeCluster(10) # running on 10 cores
      
      clusterExport(cl = c1, varlist = 'simulate_counts') # exporting functions
      
      # loading appropriate libraries
      clusterEvalQ(c1, {
        library(MASS)
      })
      
      fit <- parSapply(cl = c1,  X = samples, FUN = simulated_power,
                       params = coef_levels[ab], eff = trends[tr],
                       theta = theta_levels[ab], time = years_sampled[[ys]],
                       n_sim = 1000)
      
      stopCluster(c1)
      
      # there might be a better way to do this
      scenario_row <- which((out$Trend == names(trends[tr])) &
                              (out$Sampling == names(years_sampled[ys])) &
                              (out$Abundance == names(coef_levels[ab])))
      
      # how many samples for each power scenario?
      out$`Sample Size, Power = 0.60`[scenario_row] <-
        samples[min(which(fit >= 0.60))]
      
      out$`Sample Size, Power = 0.70`[scenario_row] <-
        samples[min(which(fit >= 0.70))]
      
      out$`Sample Size, Power = 0.80`[scenario_row] <-
        samples[min(which(fit >= 0.80))]
      
    }
  }
}

