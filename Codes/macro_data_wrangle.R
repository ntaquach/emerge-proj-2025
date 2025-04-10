# SCRIPT TO WRANGLE NEON MACRO DATA ======================

library(tidyverse)
library(ggplot2)
library(tidyr)
library(writexl)

rm(list = ls())

csv_files <- list.files("Data", pattern = "\\.csv$", full.names = TRUE)

# Prepare a list to store the combined data frames
wide_dfs <- list()
taxa_list <- list()

# Loop through each CSV file
for (file in csv_files) {
  base_name <- tools::file_path_sans_ext(basename(file))
  
  # Read and process the data
  df_wide <- read.csv(file) %>%
    select(siteID, acceptedTaxonID, estimatedTotalCount, scientificName) %>%
    group_by(siteID, acceptedTaxonID) %>%
    summarise(estimatedTotalCount = sum(estimatedTotalCount), .groups = "drop") %>%
    pivot_wider(
      names_from = acceptedTaxonID,
      values_from = estimatedTotalCount,
      values_fill = 0
    )
  
  # Assign the result to a new variable with "_wide" suffix
  assign(paste0(base_name, "_wide"), df_wide)
  
  # Extract unique taxon data and store in taxa_list
  taxa_df <- read.csv(file) %>%
    select(acceptedTaxonID, scientificName) %>%
    distinct()  # Get unique taxonID and scientificName
  
  taxa_list[[base_name]] <- taxa_df  # Store taxa info in the list
  
  # Add the wide dataframe to the wide_dfs list
  wide_dfs[[paste0(base_name, "_wide")]] <- df_wide
}

#====================MAKE A TAXA SHEET=====================

write_xlsx(list("Taxa" = bind_rows(taxa_list)), "Data/Combined NEON taxa.xlsx") #write taxa sheet for information later

#================COMBINE DATA FOR ORDINATION==========

# Get all *_wide data frames currently in your environment
wide_dfs <- mget(ls(pattern = "_wide$"))

# Combine them
combined_macro_data <- bind_rows(wide_dfs, .id = "source")

# Replace NA with 0 for taxa counts
combined_macro_data[is.na(combined_macro_data)] <- 0
combined_macro_data <- combined_macro_data %>% select(-source)

write_xlsx(combined_macro_data, "Data/Combined NEON Macro Data.xlsx")


