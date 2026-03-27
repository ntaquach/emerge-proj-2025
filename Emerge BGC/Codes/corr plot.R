##Assembling master dataset ==============================
library(ggplot2)
library(tidyverse)
library(genzplyr)
library(purrr)
library(Hmisc)
library(corrplot)
#import purrr#import Vlah data
rm(list=ls())



#load data

doc_neon <- read.csv("./Emerge BGC/Data/monthly_doc.csv") #DOC data
q_neon <- read.csv("./Emerge BGC/Data/monthly_q.csv") #discharge data
wq_neon <-  read.csv("./Emerge BGC/Data/monthly_wq.csv") #water quality data
no3_neon <- read.csv("./Emerge BGC/Data/monthly_no3.csv") #no3 data
temp_neon <- read.csv("./Emerge BGC/Data/monthly_temp.csv") #water temp data
gas_neon <- read.csv("./Emerge BGC/Data/monthly_gas.csv") #GHG data

#rename column to make joining more clear

doc_neon <- doc_neon %>% rename(doc_monthly = monthly_avg)
q_neon   <- q_neon   %>% rename(q_monthly = monthly_avg)
temp_neon<- temp_neon%>% rename(temp_monthly = monthly_avg)

wq_neon <- wq_neon %>%
  pivot_wider(names_from = var, values_from = monthly_avg) %>% 
  select(-TSS)

gas_neon <- gas_neon %>% select(c(month, year, site_code,CH4_uM, CO2_uM, N2O_uM))


#join data to make a master dataset

master_neon <- list(q_neon, wq_neon, no3_neon, temp_neon, gas_neon) %>%
  reduce(left_join, by = c("month", "year", "site_code")) %>%
  left_join(doc_neon, by = c("month", "year", "site_code")) 

write.csv(master_neon,"./Emerge BGC/Data/master_neon.csv", row.names = FALSE)

cor_data <- master_neon %>%
  select(where(is.numeric)) %>%
  select(CH4_uM, everything()) %>%
  drop_na() %>%
  cor(method = "pearson")

rc <- rcorr(as.matrix(cor_data), type = "pearson")
cor_matrix <- rc$r
p_matrix  <- rc$P

corrplot(
  cor_matrix,
  type = "lower",
  method = "color",
  col = colorRampPalette(c("#d7191c", "white", "#2c7bb6"))(200),
  tl.col = "black",
  tl.cex = 0.8,
  addCoef.col = NULL,
  p.mat = p_matrix,
  sig.level = c(0.05, 0.01, 0.001),
  insig = "label_sig",
  pch.cex = 0.8,
  pch.col = "black"
)


##RANDOM FOREST MODELING ==========================




