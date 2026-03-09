# DAW/SAW Analysis Subgroup
# Clean data file for analysis of ask 2 -- UMR
# Combined data files

# Clear environment
rm(list=ls())

# Load required packages
library(tidyverse)

# Read in compiled data file
data <- read.csv("DozerTrawl_Data.csv")

# Remove data not associated with a POI
data <- data %>% dplyr::filter(POI != "Not a POI")

# Cleaning "POI" column - rename duplicate site names
data <- data %>% mutate(POI = ifelse(POI == 
                                       "Illinois River Pools (Starved Rock Upstream to Dresden Island)", 
                                     "Illinois River Pools (Starved Rock upstream to Dresden Island)", 
                                     POI)) 

data <- data %>% mutate(POI = ifelse(POI == 
                                       "MSR pools 20-25", 
                                     "MSR Pools 20-25", 
                                     POI)) 

unique(data$POI)

################ Categorizing into abundances - medium and high ###############
# Separate data into each POI
#dataAP <- data %>% filter(POI == "Illinois River Pools (Alton upstream to Peoria)")
#dataSD <- data %>% filter(POI == "Illinois River Pools (Starved Rock upstream to Dresden Island)")
#data26 <- data %>% filter(POI == "MSR Pool 26 and below")
#data20 <- data %>% filter(POI == "MSR Pools 20-25")
#data15 <- data %>% filter(POI == "MSR Pools 15-19")
#dataDM <- data %>% filter(POI == "Des Moines River")
#dataSK <- data %>% filter(POI == "Skunk River")

# Create CPUE column in each POI dataset
#dataAP <- dataAP %>% mutate(CPUE = CountLarge / FishTime * 60)
#dataSD <- dataSD %>% mutate(CPUE = CountLarge / FishTime * 60)
#data26 <- data26 %>% mutate(CPUE = CountLarge / FishTime * 60)
#data20 <- data20 %>% mutate(CPUE = CountLarge / FishTime * 60)
#data15 <- data15 %>% mutate(CPUE = CountLarge / FishTime * 60)
#dataDM <- dataDM %>% mutate(CPUE = CountLarge / FishTime * 60)
#dataSK <- dataSK %>% mutate(CPUE = CountLarge / FishTime * 60)

# Check mean CPUE for each POI
#mean(dataAP$CPUE) # = 132.0
#mean(dataSD$CPUE) # = 34.6
#mean(data26$CPUE) # = 97.2
#mean(data20$CPUE) # = 40.9
#mean(data15$CPUE) # = 0.8
#mean(dataDM$CPUE) # = 43.1
#mean(dataSK$CPUE) # = 37.9

# Assign an abundance level to each POI
data <- data %>% 
  mutate(Abundance = case_when(
    POI %in% c("Illinois River Pools (Alton upstream to Peoria)",
               "MSR Pool 26 and below") ~
      "High",
    POI %in% c("Illinois River Pools (Starved Rock upstream to Dresden Island)",
               "MSR Pools 20-25", 
               "Des Moines River",
               "Skunk River") ~
      "Moderate",
    POI %in% c("MSR Pools 15-19") ~
      "Rare"
  ))




# Separate data into each Abundance
#dataH <- data %>% filter(Abundance == "High")
#dataM <- data %>% filter(Abundance == "Moderate")
#dataR <- data %>% filter(Abundance == "Rare")

# Create CPUE column in each Abundance dataset
#dataH <- dataH %>% mutate(CPUE = CountLarge / FishTime * 60)
#dataM <- dataM %>% mutate(CPUE = CountLarge / FishTime * 60)
#dataR <- dataR %>% mutate(CPUE = CountLarge / FishTime * 60)

# Check mean CPUE for each Abundance
#mean(dataH$CPUE) # = 120.9
#mean(dataM$CPUE) # = 36.6
#mean(dataR$CPUE) # = 0.8
