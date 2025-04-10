# SCRIPT TO WRANGLE NEON MACRO DATA ======================

library(tidyverse)
library(ggplot2)
library(tidyr)

rm(list = ls())

csv_files <- list.files("Data", pattern = "\\.csv$", full.names = TRUE)

# Loop through each CSV file
for (file in csv_files) {
  # Get the same name
  base_name <- tools::file_path_sans_ext(basename(file))
  
  # Read and process the data
  df_wide <- read.csv(file) %>%
    select(siteID, acceptedTaxonID, estimatedTotalCount) %>%
    group_by(siteID, acceptedTaxonID) %>%
    summarise(estimatedTotalCount = sum(estimatedTotalCount), .groups = "drop") %>%
    pivot_wider(
      names_from = acceptedTaxonID,
      values_from = estimatedTotalCount)
  # Assign the result to a new variable with "_wide" suffix
  assign(paste0(base_name, "_wide"), df_wide)
}

rm(df_wide)

str(ARIK_SPRING_wide)

# Get all *_wide data frames currently in your environment
wide_dfs <- mget(ls(pattern = "_wide$"))

# Combine them
combined_df <- bind_rows(wide_dfs, .id = "source")

# Replace NA with 0 for taxa counts
combined_df[is.na(combined_df)] <- 0
combined_df <- combined_df %>% select(-source)





