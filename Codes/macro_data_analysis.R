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

neon_macro_only <- neon_macro[,2:ncol(neon_macro)] #exclude the first columns that contains SU information


skewness(neon_macro_only)
kurtosis(neon_macro_only) 

#NMDS is better because skewness and kurtosis of data are too far from 0.

#=============RUN NMDS============

set.seed(1)
neon_macro_nmds <- metaMDS(neon_macro_only,distance="bray",maxit=999,trymax = 500,wascores = T,k=2,autotransform = F)

#Plot NMDS

data.scores <- as.data.frame(scores(neon_macro_nmds)$sites) #extract NMDS scores
data.scores$siteID <- neon_macro$siteID

neon_macro_graph <- ggplot(data.scores, aes(NMDS1, NMDS2)) + 
  geom_point(aes(shape = siteID),size = 5)+ 
  scale_shape_manual(values = 1:10) + # Manually set shapes (1 to 10)
  labs(
    title = "",
    x = "NMDS1", shape = "NEON Sites", y = "NMDS2")  + 
  geom_hline(yintercept=0, linetype="dashed", 
             color = "black", size=1) + #add horizontal and vertical lines at 0
  geom_vline(xintercept=0, linetype="dashed", 
             color = "black", size=1) + # +geom_text(aes(label=Site),hjust=-0.25, vjust=1.5,size=4); option for sites
  anh_theme()

neon_macro_graph

ggsave("NEON Macro Spring NMDS.png",neon_macro_graph,width=8,height=5)
