library(tidyverse)

# Read in compiled data file
UMRdata_clean <- read.csv("UMRdata2.csv")

# Remove data not associated with a POI
UMRdata_clean <- UMRdata_clean %>% 
  filter(POI != "Not a POI")

# Standardize POI names
UMRdata_clean <- UMRdata_clean %>%
  mutate(
    POI = ifelse(
      POI == "Illinois River Pools (Starved Rock Upstream to Dresden Island)",
      "Illinois River Pools (Starved Rock upstream to Dresden Island)",
      POI
    ),
    POI = ifelse(
      POI == "MSR pools 20-25",
      "MSR Pools 20-25",
      POI
    )
  )

# Assign abundance categories
UMRdata_clean <- UMRdata_clean %>%
  mutate(Abundance = case_when(
    POI %in% c(
      "Illinois River Pools (Alton upstream to Peoria)",
      "MSR Pool 26 and below"
    ) ~ "High",
    POI %in% c(
      "Illinois River Pools (Starved Rock upstream to Dresden Island)",
      "MSR Pools 20-25",
      "Des Moines River",
      "Skunk River"
    ) ~ "Moderate",
    POI %in% c("MSR Pools 15-19") ~ "Rare",
    TRUE ~ NA_character_
  ))

# Remove samples with zero/negative FishTime
UMRdata_clean <- UMRdata_clean %>% 
  filter(FishTime > 0)

