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

neon_macro <- read_xlsx("Data/Combined NEON Macro Data.xlsx")
site_info <- read_xlsx("Data/Site Type.xlsx")

neon_macro_site <- site_info %>% left_join(neon_macro, by = "siteID")
neon_macro_only <- neon_macro[,4:ncol(neon_macro)] #exclude the first columns that contains SU information


skewness(neon_macro_only)
kurtosis(neon_macro_only) 

#NMDS is better because skewness and kurtosis of data are too far from 0.

#=============RUN NMDS============

set.seed(1)
neon_macro_nmds <- metaMDS(neon_macro_only,distance="bray",maxit=999,trymax = 500,wascores = T,k=2,autotransform = F)

#Plot NMDS

data.scores <- as.data.frame(scores(neon_macro_nmds)$sites) #extract NMDS scores
data.scores$siteID <- neon_macro_site$siteID
data.scores$siteType <- as.factor(neon_macro_site$siteType)
data.scores$domain <- as.factor(neon_macro_site$domain)

#====set colors===
domain_colors <- c(
  "#E41A1C", "#377EB8", "#4DAF4A", "#FF7F00", "#FFFF33", 
  "#A65628", "#F781BF", "#999999", "#66C2A5", "#FC8D62", 
  "#8DA0CB", "#E78AC3"  # 12 distinct colors
)

#plot  graph==============

neon_macro_graph <- ggplot(data.scores, aes(NMDS1, NMDS2)) + 
  geom_point(aes(shape = siteType, fill=domain),size = 5, color="black")+ 
  scale_shape_manual(values = c(21,23)) + # Manually set shapes
  scale_fill_manual(values = domain_colors)+
  labs(
    title = "",
    x = "NMDS1", shape = "NEON Sites", y = "NMDS2", fill="Domain")  + 
  geom_hline(yintercept=0, linetype="dashed", 
             color = "black", size=1) + #add horizontal and vertical lines at 0
  geom_vline(xintercept=0, linetype="dashed", 
             color = "black", size=1) + # +geom_text(aes(label=Site),hjust=-0.25, vjust=1.5,size=4); option for sites
  guides(fill = guide_legend(override.aes = list(shape = 21, color = "black"))) +
  geom_text(aes(label=siteID),hjust=0.4, vjust=2,size=3, color="black") +
  anh_theme()

neon_macro_graph

ggsave("Figures/NEON Macro Spring NMDS.png",neon_macro_graph,width=8,height=5)

#=========extract species cor with NMDS axis

species_axis_cor <- cor(neon_macro_only, data.scores, method = "pearson")
species_cor_df <- data.frame(Species = rownames(species_axis_cor),
                             species_axis_cor,
                             row.names = NULL)
write_xlsx(species_cor_df, "species_nmds_correlations.xlsx")

#======extract spp score

taxon.scores <- as.data.frame(scores(neon_macro_nmds,display=c("species")))
taxa_nmds1 <- taxon.scores[c('GLOSP1','SPESP4','HYDSP26','ABLSP','NAISP','NAICOM','CHISP6'),]
taxa_nmds2 <- taxon.scores[c('CHISP8','PRIPRO','HYDSP21','HETCOR','HYDSP18'),]

#plot biplot=============
neon_biplot <- neon_macro_graph + 
  # Add stars for sensitive taxa
  geom_point(data = taxa_nmds1, aes(x = NMDS1, y = NMDS2), shape = 8, size = 2, color = "black") +
  geom_text(data = taxa_nmds1,
            aes(x = NMDS1, y = NMDS2, label = rownames(taxa_nmds1),
                hjust = 0.25 * (1 - sign(NMDS1)), vjust = 0.5 * (1 - sign(NMDS2))),
            color = "black", size = 4, fontface = "bold") +
  
  # Add stars for tolerant taxa
  geom_point(data = taxa_nmds2, aes(x = NMDS1, y = NMDS2), shape = 9, size = 2, color = "black") +
  geom_text(data = taxa_nmds2,
            aes(x = NMDS1, y = NMDS2, label = rownames(taxa_nmds2),
                hjust = 0.25 * (1 - sign(NMDS1)), vjust = 0.5 * (1 - sign(NMDS2))),
            color = "black", size = 4, fontface = "italic")

neon_biplot
#=========
ggsave("Figures/NEON NMDS Biplot.png",neon_biplot,width=9,height=5)
