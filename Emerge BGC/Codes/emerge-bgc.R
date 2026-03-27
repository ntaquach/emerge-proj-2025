library(ggplot2)
library(tidyverse)
library(genzplyr)
#import Vlah data
rm(list=ls())
neon<-read.csv("./Emerge BGC/Data/timeseries_neon.csv") #load neon data

unique(neon$var) #see what variables are available

neon$date <- as.Date(neon$date) #make sure the format of the date is Date

                              

                                  #working with DOC data ===================================

doc_neon <- neon %>% filter(date > "2021-01-01" & date < "2022-12-31" & var == "DOC"
                            & ms_status == 1 & ms_interp == 1 & var_category == "stream_chemistry") #export DOC data
# 
# # you can change to whichever variable you are interested in by changing var and var category
# #ms_status and ms_interp are QA QC and I believe 0 represents good data 
# # 
# unique(doc_neon$site_code) #check out what sites are available
# 
# write.csv(doc_neon,"./Emerge BGC/Data/doc.csv", row.names = FALSE) #write our csv dataset

doc_neon <- read.csv("./Emerge BGC/Data/doc.csv")

ggplot(doc_neon, aes(x=date, y=val)) + geom_point() +
  facet_wrap(~site_code, scales="free_y")

  #============summarize DOC data
    #====monthly=======

  monthly_doc_neon <- doc_neon %>% 
    glow_up(month = lubridate::month(mdy(date)),
            year = lubridate::year(mdy(date))) %>%
    squad_up(month, year, site_code) %>%
    no_cap(monthly_avg = mean(val))

write.csv(monthly_doc_neon,"./Emerge BGC/Data/monthly_doc.csv", row.names = FALSE)

ggplot(monthly_doc_neon, aes(x=factor(year), y=monthly_avg)) + geom_boxplot() +
  facet_wrap(~site_code, scales="free_y")
    
    #seasonally=========
  season_doc_neon <- doc_neon %>% 
    glow_up(month = lubridate::month(mdy(date)),
            year = lubridate::year(mdy(date))) %>%
    mutate(
      season = case_when(
        month %in% c(12, 1, 2) ~ "Winter",
        month %in%  c(3,4,5)  ~ "Spring",
        month %in%  c(6,7,8)  ~ "Summer",
        month %in%  c(9,10,11)  ~ "Fall")) %>%
    squad_up(season, year, site_code) %>%
    no_cap(season_avg = mean(val))

write.csv(season_doc_neon,"./Emerge BGC/Data/season_doc.csv", row.names = FALSE)


ggplot(season_doc_neon, aes(x=factor(year), y=season_avg, color=season)) + 
  geom_point() +
  facet_wrap(~site_code, scales="free_y")


                                              #=====Discharge data===============

Q_neon <- neon %>% filter(date > "2021-01-01" & date < "2022-12-31" & var == "discharge"
                            & ms_status == 0 & ms_interp == 0 & var_category == "discharge")
# unique(Q_neon$site_code)
write.csv(Q_neon,"./Emerge BGC/Data/Q.csv", row.names = FALSE)

q_neon <- read.csv("./Emerge BGC/Data/Q.csv")

ggplot(q_neon, aes(x=date, y=val)) + geom_point() +
  facet_wrap(~site_code, scales="free_y")


  #===summarize Q data

#====monthly=======

monthly_q_neon <- q_neon %>% 
  glow_up(month = lubridate::month(ymd(date)),
          year = lubridate::year(ymd(date))) %>%
  squad_up(month, year, site_code) %>%
  no_cap(monthly_avg = mean(val))

write.csv(monthly_q_neon,"./Emerge BGC/Data/monthly_q.csv", row.names = FALSE)

ggplot(monthly_q_neon, aes(x=factor(year), y=monthly_avg)) + geom_boxplot() +
  facet_wrap(~site_code, scales="free_y")

#seasonally=========
season_q_neon <- q_neon %>% 
  glow_up(month = lubridate::month(ymd(date)),
          year = lubridate::year(ymd(date))) %>%
  mutate(
    season = case_when(
      month %in% c(12, 1, 2) ~ "Winter",
      month %in%  c(3,4,5)  ~ "Spring",
      month %in%  c(6,7,8)  ~ "Summer",
      month %in%  c(9,10,11)  ~ "Fall")) %>%
  squad_up(season, year, site_code) %>%
  no_cap(season_avg = mean(val))

write.csv(season_q_neon,"./Emerge BGC/Data/season_q.csv", row.names = FALSE)


ggplot(season_q_neon, aes(x=factor(year), y=season_avg, color=season)) + 
  geom_point() +
  facet_wrap(~site_code, scales="free_y")


#===summarize water quality data

# wq_neon <- neon %>% filter(date > "2021-01-01" & date < "2022-12-31" & var == c("pH", "spCond", "DO",
#                                                                                 "DO_sat", "PAR", "TSS")
#                             & ms_status == 0 & ms_interp == 0 & var_category == "stream_chemistry") #export water quality data
# 
write.csv(wq_neon,"./Emerge BGC/Data/wq.csv", row.names = FALSE) #write our csv dataset

wq_neon <- read.csv("./Emerge BGC/Data/wq.csv")
#====monthly=======

monthly_wq_neon <- wq_neon %>% 
  glow_up(month = lubridate::month(ymd(date)),
          year = lubridate::year(ymd(date))) %>%
  squad_up(month, year, site_code, var) %>%
  no_cap(monthly_avg = mean(val))

write.csv(monthly_wq_neon,"./Emerge BGC/Data/monthly_wq.csv", row.names = FALSE)

ggplot(monthly_wq_neon, aes(x=factor(year), y=monthly_avg)) + geom_boxplot() +
  facet_wrap(~interaction(site_code,var), scales="free_y")

#seasonally=========
season_wq_neon <- wq_neon %>% 
  glow_up(month = lubridate::month(ymd(date)),
          year = lubridate::year(ymd(date))) %>%
  mutate(
    season = case_when(
      month %in% c(12, 1, 2) ~ "Winter",
      month %in%  c(3,4,5)  ~ "Spring",
      month %in%  c(6,7,8)  ~ "Summer",
      month %in%  c(9,10,11)  ~ "Fall")) %>%
  squad_up(season, year, site_code,var) %>%
  no_cap(season_avg = mean(val))

write.csv(season_wq_neon,"./Emerge BGC/Data/season_wq.csv", row.names = FALSE)


ggplot(season_wq_neon, aes(x=factor(year), y=season_avg, color=season)) + 
  geom_point() +
  facet_wrap(~interaction(site_code,var), scales="free_y")
