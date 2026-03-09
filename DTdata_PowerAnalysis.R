# DAW/SAW Analysis Subgroup
# Power Analysis for Upper Mississippi River Sub-basin
# Sampling effort needed for:
#     60, 70, & 80% power 
#     to detect a - 15, 25, & 50% change in 
#     adult Silver Carp abundance over 
#     a 3- and 5-year time span at 
#     alpha = 0.05

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


######################## Define Functions #####################################

# Simulate count data from NB model
# applies offset & effect size to get new expected count
# returns data frame with counts over time

simulate_counts <- function(n, params, eff, theta, offset = log(5),
                            time = 1:3) {
  # n - sample size
  # params - coefficients of nb_all model
  # eff - effect size (influence of params on outcome)
  # theta - dispersion from model
  # offset - effort offset (log FishTime), multiplies by standard sampling
    # defaults to 5 minutes for dozer trawl data
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

simulated_power <- function(sample_size, n_sim = 10000, params, eff, theta,
                            alpha = 0.05, offset = log(5), time){
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
samples <- seq(25, 500, by = 5)

# defining simulation scenarios
trends <- c(log(0.85), log(0.75), log(0.50))
names(trends) <- c('15% Reduction', '25% Reduction', '50% Reduction')

years_sampled <- list(1:3, 1:5, c(1, 3, 5))
names(years_sampled) <- c('Annually for 3 years', 'Annually for 5 years',
                          'Year 1, 3, and 5')

coef_levels <- c(coef(nb_high), coef(nb_mod))#, coef(nb_rare))
names(coef_levels) <- c('High Abundance', 'Moderate Abundance')#, 'Rare Abundance')

theta_levels <- c(nb_high$theta, nb_mod$theta)#, nb_rare$theta)
names(theta_levels) <- c('High Abundance', 'Moderate Abundance')#, 'Rare Abundance')

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
                       n_sim = 100)
      
      stopCluster(c1)
      
      # there might be a better way to do this
      scenario_row <- which((out$Trend == names(trends[tr])) &
                              (out$Sampling == names(years_sampled[ys])) &
                              (out$Abundance == names(coef_levels[ab])))
      
      # how many samples for each power scenario?
      #out$`Sample Size, Power = 0.60`[scenario_row] <-
      #  samples[min(which(fit >= 0.60))]
      
      #out$`Sample Size, Power = 0.70`[scenario_row] <-
      #  samples[min(which(fit >= 0.70))]
      
      #out$`Sample Size, Power = 0.80`[scenario_row] <-
      #  samples[min(which(fit >= 0.80))]
      
      out$`Sample Size, Power = 0.60`[scenario_row] <-
        if (length(which(fit >= 0.60)) > 0) {
          samples[min(which(fit >= 0.60))]
        } else {
          NA
        }
      
      out$`Sample Size, Power = 0.70`[scenario_row] <-
        if (length(which(fit >= 0.70)) > 0) {
          samples[min(which(fit >= 0.70))]
        } else {
          NA
        }
      
      out$`Sample Size, Power = 0.80`[scenario_row] <-
        if (length(which(fit >= 0.80)) > 0) {
          samples[min(which(fit >= 0.80))]
        } else {
          NA
        }
      
    }
  }
}





## THIS IS WHERE ROTA STOPPED MUCKING AROUNG


simulate_power <- function(n_sim, sample_size, time, params, theta, eff,
                           alpha = 0.05, offset = log(5), num_cores = 3) {

  
  c1 <- makeCluster(num_cores)
  
  
  
  
  
  power_results <- parLapply(c1, abun, function(abundance) {
    library(dplyr)
    
    results_list <- data.frame(sample_size = numeric(0),
                               power_DT = numeric(0),
                               Abundance = character())

          for (size in sample_size) {
          p_val_DT <- numeric(n_sim)
          
          for (i in 1:n_sim) {
            count_DT <- simulate_counts(n = size,
                                        params = params,
                                        eff = eff,
                                        theta = theta,
                                        offset = offset,
                                        time = time,
                                        abun = abundance)
            
            # testing for warnings or errors from the fitted model
            test_DT <- try(assertCondition(glm.nb(y ~ Year + offset(log_offset),
                                                  data = count_DT)), silent = T)
            
            # simulating data until no warnings
            while(any(class(test_DT[[1]]) %in%
                      c('simpleWarning', 'simpleError'))){
              count_DT <- simulate_counts(n = size,
                                          params = params,
                                          eff = eff,
                                          theta = theta,
                                          offset = offset,
                                          time = time,
                                          abun = abundance)
              
              # testing for warnings or errors from the fitted model
              test_DT <-
                try(assertCondition(glm.nb(y ~ Year + offset(log_offset),
                                           data = count_DT)), silent = T)
            }
            
            modDT <- glm.nb(y ~ Year + offset(log_offset), data = count_DT)
            
            p_val_DT[i] <- coef(summary(modDT))["Year", "Pr(>|z|)"]
            
          }
          
          results_list <- bind_rows(results_list, 
                                    data.frame(Abundance = abundance,
                                               sample_size = size,
                                               power_DT = mean(p_val_DT < 0.05)))
          }
    
    return(results_list)
  })
  
  stopCluster(c1)
  return(bind_rows(power_results))
}


########################## Run power simulation ###############################
# Adjust effect size: 15% decrease in abundance over 3 year
eff <- log(0.85)/3

y3_0.85 <- simulate_power(n_sim = 1000, 
                       #sample_size = seq(1750, 2150, by = 50), # 60
                       #sample_size = seq(2400, 2750, by = 50), # 70
                       sample_size = seq(3250, 3350, by = 50), # 80
                       #sample_size = seq(1750, 3350, by = 50), # combined
                       time = 1:3, 
                       params = params, 
                       theta = theta, 
                       eff = eff, 
                       abun = #c(
                         'Moderate')#, 'High')) # Moderate or High

saveRDS(y3_0.85, '3y-0.85.rds')
y3_0.85 <- readRDS('3y-0.85.rds')

y3_0.85_60p <- y3_0.85 %>%
  filter(power_DT >= 0.6) %>%
  group_by(Abundance) %>%
  summarize(MinSites = min(sample_size))

print(y3_0.85_60p)

y3_0.85_70p <- y3_0.85 %>%
  filter(power_DT >= 0.7) %>%
  group_by(Abundance) %>%
  summarize(MinSites = min(sample_size))

print(y3_0.85_70p)

y3_0.85_80p <- y3_0.85 %>%
  filter(power_DT >= 0.8) %>%
  group_by(Abundance) %>%
  summarize(MinSites = min(sample_size))

print(y3_0.85_80p)


# Adjust effect size: 25% decrease in abundance over 3 year
eff <- log(0.75)/3

y3_0.75 <- simulate_power(n_sim = 1000, 
                          sample_size = seq(750, 2500, by = 50),
                          time = 1:3,
                          params = params,
                          theta = theta,
                          eff = eff)

saveRDS(y3_0.75, '3y-0.75.rds')
y3_0.75 <- readRDS('3y-0.75.rds')

y3_0.75_60p <- y3_0.75 %>%
  filter(power_DT >= 0.6) %>%
  group_by(Abundance) %>%  
  summarize(MinSites = min(sample_size))

print(y3_0.75_60p)

y3_0.75_70p <- y3_0.75 %>%
  filter(power_DT >= 0.7) %>%
  group_by(Abundance) %>%  
  summarize(MinSites = min(sample_size))

print(y3_0.75_70p)

y3_0.75_80p <- y3_0.75 %>%
  filter(power_DT >= 0.8) %>%
  group_by(Abundance) %>%  
  summarize(MinSites = min(sample_size))

print(y3_0.75_80p)


# Adjust effect size: 50% decrease in abundance over 3 year
eff <- log(0.5)/3

y3_0.5 <- simulate_power(n_sim = 1000, 
                         sample_size = seq(750, 2500, by = 50),
                         time = 1:3,
                         params = params,
                         theta = theta,
                         eff = eff)

saveRDS(y3_0.5, '3y-0.5.rds')
y3_0.5 <- readRDS('3y-0.5.rds')

y3_0.5_60p <- y3_0.5[[1]] %>%
  filter(power_DT >= 0.6) %>%
  group_by(Abundance) %>%  
  summarize(MinSites = min(sample_size))

print(y3_0.5_60p)

y3_0.5_70p <- y3_0.5[[1]] %>%
  filter(power_DT >= 0.7) %>%
  group_by(Abundance) %>%  
  summarize(MinSites = min(sample_size))

print(y3_0.5_70p)

y3_0.5_80p <- y3_0.5[[1]] %>%
  filter(power_DT >= 0.8) %>%
  group_by(Abundance) %>%  
  summarize(MinSites = min(sample_size))

print(y3_0.5_80p)
