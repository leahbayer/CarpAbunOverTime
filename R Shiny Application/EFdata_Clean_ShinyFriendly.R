library(tidyverse)

# Read in compiled data file
ILRdataEF_clean <- read.csv("ILRdataEF.csv")

# Remove data not associated with a POI
ILRdataEF_clean <- ILRdataEF_clean %>% 
  filter(POI != "Not a POI")

# Standardize POI name
ILRdataEF_clean <- ILRdataEF_clean %>%
  mutate(POI = ifelse(
    POI == "Illinois River Pools (Starved Rock Upstream to Dresden Island)",
    "Illinois River Pools (Starved Rock upstream to Dresden Island)",
    POI
  ))

# Remove low-abundance regions
ILRdataEF_clean <- ILRdataEF_clean %>%
  filter(!POI %in% c("MSR Pools 15-19", "MSR Pool 5a-8"))

# Assign abundance category
ILRdataEF_clean <- ILRdataEF_clean %>%
  mutate(Abundance = "Moderate")

# Remove samples with zero/negative FishTime
ILRdataEF_clean <- ILRdataEF_clean %>% 
  filter(FishTime > 0)