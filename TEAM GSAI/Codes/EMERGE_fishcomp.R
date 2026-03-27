# Load necessary libraries
library(dplyr)
library(tidyr)
library(lubridate)
library(vegan)
library(ggplot2)
library(readxl)
library(moments) #load moments package to check for kurtosis and skewness

file_paths <- list.files(path = "~/Downloads", pattern = "fsh_perFish.*\\.csv$", full.names = TRUE)
file_paths1 <- list.files(path = "~/Downloads", pattern = "fsh_fieldData.*\\.csv$", full.names = TRUE)

# Define columns to extract from fish data
fish_columns <- c("siteID","taxonID", "scientificName", "domainID", "boutID")

# Define columns to extract from field data
field_columns <- c("boutID", "domainID", "siteID")
# Function to load and select relevant columns from fish data
load_fish_data <- function(file) {
  read.csv(file) %>%
    select(any_of(fish_columns))
}

# Function to load and select relevant columns from field data
load_field_data <- function(file) {
  read.csv(file) %>%
    select(any_of(field_columns))
}

# Read and combine fish data
fish_df <- file_paths %>%
  lapply(load_fish_data) %>%
  bind_rows()

# Read and combine field data
field_df <- file_paths1 %>%
  lapply(load_field_data) %>%
  bind_rows()

# Merge fish and field data
combined_df <- fish_df %>%
  left_join(field_df, by = "domainID")

# Create a wide matrix with species as columns and count per site and season
fish_all <- combined_df %>%
  group_by(siteID.x, boutID, scientificName, domainID) %>%
  summarise(count = n(), .groups = "drop") %>%
  pivot_wider(names_from = scientificName, values_from = count, values_fill = 0)
names(fish_all) <- gsub(" ", "_", names(fish_all))
names(fish_all)[names(fish_all)=="siteID.x"] <- "siteID"
names(fish_all)[names(fish_all)=="boutID"] <- "season"

# Count how many sites have data per season
site_season_summary <- fish_all %>%
  filter(!is.na(season)) %>%
  distinct(siteID, season) %>%
  count(season, name = "n_sites") %>%
  arrange(desc(n_sites))
# List of each site and their respective seasons
site_season_list <- fish_all %>%
  filter(!is.na(season)) %>%
  distinct(siteID, season) %>%
  arrange(siteID, season)

#
abund_spr <- fish_all %>% 
  filter(season == "spring") 

abund_mat <- abund_spr %>%
  select(-siteID, -season,-domainID)

skewness(abund_mat)
kurtosis(abund_mat) 

# 4. NMDS on abundance (Bray)
set.seed(1)
neon_fish_nmds <- metaMDS(abund_mat, distance="bray",maxit=999,trymax = 500,wascores = T,k=2,autotransform = F)
#Plot NMDS
data.scores <- as.data.frame(scores(neon_fish_nmds)$sites) #extract NMDS scores
data.scores$siteID <- abund_spr$siteID
data.scores$domain <- as.factor(abund_spr$domainID)

#plot  graph==============

neon_fish_graph <- ggplot(data.scores, aes(NMDS1, NMDS2,fill = domain)) + 
  geom_point(shape=21,size = 5,alpha = 0.8)+ 
  labs(
    title = "",
    x = "NMDS1", shape = "NEON Sites", y = "NMDS2", fill="Domain")  + 
  geom_hline(yintercept=0, linetype="dashed", 
             color = "black", size=1) + #add horizontal and vertical lines at 0
  geom_vline(xintercept=0, linetype="dashed", 
             color = "black", size=1) + 
  guides(fill = guide_legend(override.aes = list(shape = 21, color = "black"))) +
  geom_text(aes(label=siteID),hjust=0.4, vjust=2,size=3, color="black") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))


# Save NMDS plot
ggsave("NMDS_fish_abundance.png", neon_fish_graph, width = 6, height = 5, dpi = 600)

#####Environmental data
# Load environmental data
env_df <- read.csv("~/Downloads/env_seasonal.csv")%>%
  filter(season == "Spring")
env_df$season <- tolower(env_df$season)

# Final combined dataset: environment and fish data
final_combined <- inner_join(env_df, fish_all, by = c("siteID", "season")) %>%
  select(-("X"), -("domainID.x"),-("domainID.y"))

# Save result 
write.csv(final_combined, "final_combined_fish_env.csv", row.names = FALSE)


# Perform PCoA on abundance data
abundance_data <- fish_all %>%
  select(-siteID, -season,-domainID)

# Compute Bray-Curtis dissimilarity and run PCoA
bray_dist <- vegdist(abundance_data, method = "bray")
pcoa_result <- cmdscale(bray_dist, eig = TRUE, k = 2)

# Calculate percent variance explained by axes
var_explained <- round(100 * pcoa_result$eig / sum(pcoa_result$eig), 1)

# Create a PCoA dataframe for ggplot
pcoa_df <- data.frame(SiteID = fish_all$siteID,
                      Season = fish_all$season,
                      Domain = fish_all$domainID,
                      PC1 = pcoa_result$points[, 1],
                      PC2 = pcoa_result$points[, 2])

# Save 
write.csv(pcoa_df, "scores_pcoa_fish.csv", row.names = FALSE)

# Plot PCoA using ggplot2 with explained variance for Season
ggplot(pcoa_df, aes(x = PC1, y = PC2, color = Season)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_minimal() +
  labs(title = "PCoA of Fish Community Composition",
       x = paste0("PCoA Axis 1 (", var_explained[1], "%)"),
       y = paste0("PCoA Axis 2 (", var_explained[2], "%)"),
       color = "Season") +
  theme(plot.title = element_text(hjust = 0.5))

# Plot PCoA using ggplot2 with explained variance for Domain
pcoa_plot<-ggplot(pcoa_df, aes(x = PC1, y = PC2, color = Domain, shape = Season)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_minimal() +
  geom_text(aes(label = SiteID), size = 3, vjust = -0.5, show.legend = FALSE) +
  labs(title = "PCoA of Fish Community Composition",
       x = paste0("PCoA Axis 1 (", var_explained[1], "%)"),
       y = paste0("PCoA Axis 2 (", var_explained[2], "%)"),
       color = "DomainID") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

ggsave("PCoA_fish_composition.png", plot = pcoa_plot, dpi = 600, width = 8, height = 6)