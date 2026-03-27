##================SCRIPT TO VISUALIZE NEON MACRO COMMUNITY COMPOSITION===============

library(readxl)
library(vegan)
library(ggplot2)
library(tidyverse)
library(moments) #load moments package to check for kurtosis and skewness


rm(list = ls())

#==============THEME FUNCTION
anh_theme <- function() {
  theme(
    axis.text.x = element_text(hjust = 0.5, colour = "black", face = "bold", size = 8), 
    axis.text.y = element_text(colour = "black", size = 8, face = "bold"), 
    legend.text = element_text(size = 8, face = "bold", colour = "black"), 
    legend.position = "right", 
    axis.title.y = element_text(face = "bold", size = 7), 
    axis.title.x = element_text(face = "bold", size = 9, colour = "black", margin = margin(t = 5)), 
    axis.title.y.left = element_text(face = "bold", size = 9, colour = "black", margin = margin(r = 10)),
    legend.title = element_text(size = 9, colour = "black", face = "bold"), 
    panel.background = element_blank(), 
    panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
    legend.key = element_blank(),
    plot.title = element_text(color = "black", size = 9, face = "bold", hjust = 0.5),
    strip.text = element_text(size = 9)
  )
}


#Data
neon_macro <- read_xlsx("Data/Combined NEON Macro Data.xlsx")
site_info <- read_xlsx("Data/Site Type.xlsx")

neon_macro_site <- site_info %>% left_join(neon_macro, by = "siteID")
neon_macro_only <- neon_macro[,4:ncol(neon_macro)] #exclude the first columns that contains SU information

# Step 2: Hellinger transformation (recommended due to dominance)
data_hel <- decostand(neon_macro_only, method = "hellinger")

# Step 3: Run PCoA on Euclidean distance (works well post-Hellinger)
dist_hel <- vegdist(data_hel, method = "euclidean")
pcoa_hel <- cmdscale(dist_hel, eig = TRUE, k = 2)

##PERMANOVA of macro data:

adonis2(dist_hel ~ domain, data = neon_macro_site,
        permutations = 999, by = "margin")

# Step 4: Extract axis values and percent variance explained
points_df <- as.data.frame(pcoa_hel$points)
colnames(points_df) <- c("PCoA1", "PCoA2")
points_df$Site <- rownames(points_df)

# Percent variance explained
eig_vals <- pcoa_hel$eig
var_explained <- eig_vals / sum(eig_vals)
pc1_var <- round(var_explained[1] * 100, 1)
pc2_var <- round(var_explained[2] * 100, 1)

# Step 5: Plot in ggplot2
#add site details back to points_Df

points_df$siteID <- neon_macro$siteID
points_df$siteType <- as.factor(site_info$siteType)
points_df$domain <- as.factor(site_info$domain)

site_plot <- ggplot(points_df, aes(x = PCoA1, y = PCoA2, label = siteID)) +
  geom_point(aes(fill=domain),size = 5, shape = 21)+ 
  geom_text(vjust = -0.5, hjust = 0.6) +
  anh_theme() +
  labs(
    title = "PCoA of Community Composition (Hellinger + Euclidean)",
    x = paste0("PCoA1 (", pc1_var, "%)"),
    y = paste0("PCoA2 (", pc2_var, "%)"),
    fill="Domain"
  ) +
  geom_hline(yintercept=0, linetype="dashed", 
             color = "black", size=0.5) + #add horizontal and vertical lines at 0
  geom_vline(xintercept=0, linetype="dashed", 
             color = "black", size=0.5) # +geom_text(aes(label=Site),hjust=-0.25, vjust=1.5,size=4); option for sites
site_plot

points_df$siteID <- neon_macro_site$siteID  # make sure this matches row order

#Using cor with axes=======================
# Calculate species-PCoA axis correlations
species_cor <- cor(data_hel, pcoa_hel$points)  # Correlate species x PCoA axes
species_cor <- as.data.frame(species_cor)      # Convert to dataframe
colnames(species_cor) <- c("PCoA1", "PCoA2")   # Rename columns
species_cor$Species <- rownames(species_cor)   # Add species names

filtered_species <- species_cor[abs(species_cor$PCoA1) > 0.5 &
                                  abs(species_cor$PCoA2) > 0.31, ]

#BIPLOT===========
pcoa_biplot <- site_plot +
  # Add species points (red)
  geom_point(
    data = filtered_species,
    aes(x = PCoA1, y = PCoA2),
    color = "red",
    size = 2,
    inherit.aes = FALSE  # Ignore main plot aesthetics
  ) +
  # Add species labels (red)
  geom_text(
    data = filtered_species,
    aes(x = PCoA1, y = PCoA2, label = Species),
    color = "black",
    vjust = 1.5,
    size = 3,
    inherit.aes = FALSE
  ) 
pcoa_biplot

ggsave("Figures/PCOA Biplot.svg", pcoa_biplot, width=6, height=6)
