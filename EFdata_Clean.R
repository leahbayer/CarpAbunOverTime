# DAW/SAW Analysis Subgroup
# Clean data file for analysis of ask 2 -- ILR EF data

# Clear environment
rm(list=ls())

# Load required packages
library(tidyverse)

# Read in compiled data file
data <- read.csv("Electrofishing_Data.csv")

# Remove data not associated with a POI
data <- data %>% dplyr::filter(POI != "Not a POI")

# Cleaning "POI" column - rename duplicate site names
data <- data %>% mutate(POI = ifelse(POI == 
                                       "Illinois River Pools (Starved Rock Upstream to Dresden Island)", 
                                     "Illinois River Pools (Starved Rock upstream to Dresden Island)", 
                                     POI)) 
#unique(data$POI)
################ Categorizing into abundances - medium and high ###############
# Separate data into each POI
#dataAP <- data %>% filter(POI == "Illinois River Pools (Alton upstream to Peoria)")
#dataSD <- data %>% filter(POI == "Illinois River Pools (Starved Rock upstream to Dresden Island)")
#data26 <- data %>% filter(POI == "MSR Pool 26 and below")
#data20 <- data %>% filter(POI == "MSR pools 20-25")
#data15 <- data %>% filter(POI == "MSR Pools 15-19")
#data5 <- data %>% filter(POI == "MSR Pool 5a-8")

# Create CPUE column in each POI dataset
#dataAP <- dataAP %>% mutate(CPUE = CountLarge / FishTime * 60)
#dataSD <- dataSD %>% mutate(CPUE = CountLarge / FishTime * 60)
#data26 <- data26 %>% mutate(CPUE = CountLarge / FishTime * 60)
#data20 <- data20 %>% mutate(CPUE = CountLarge / FishTime * 60)
#data15 <- data15 %>% mutate(CPUE = CountLarge / FishTime * 60)
#data5 <- data5 %>% mutate(CPUE = CountLarge / FishTime * 60)

# Check mean CPUE for each POI
#mean(dataAP$CPUE) # = 14.1
#mean(dataSD$CPUE) # = 6.1
#mean(data26$CPUE) # = 4.9
#mean(data20$CPUE) # = 4.8
#mean(data15$CPUE) # = 0
#mean(data5$CPUE) # = 0

data <- data %>%
  dplyr::filter(POI != "MSR Pools 15-19") %>%
  dplyr::filter(POI != "MSR Pool 5a-8")

# Assign an abundance level to each POI
# Low = < 2
# Median = 2-40
# High = > 40
data <- data %>% 
  mutate(Abundance = case_when(
    POI %in% c("Illinois River Pools (Alton upstream to Peoria)",
               "Illinois River Pools (Starved Rock upstream to Dresden Island)",
               "MSR Pool 26 and below", 
               "MSR pools 20-25") ~
      "Median"
  ))

# Remove samples with negative or zero sampling time
data <- data %>% filter(FishTime > 0)