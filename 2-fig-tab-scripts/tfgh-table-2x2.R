#######################################
# WASH Benefits Bangladesh STH KK qPCR validation
# 2x2 table of concordance discordance

# by Jade Benjamin-Chung
# jadebc@berkeley.edu
#######################################
rm(list=ls())

# configure directories, load libraries and base functions
source(paste0(here::here(), "/0-config.R"))

# load data
load(paste0(data_dir,"qdata.RData"))
load(paste0(data_dir,"kappa_test.RData"))

#---------------------------------------
# preprocess data
#---------------------------------------
qdata = qdata %>%
  mutate(positive.Al2.lab = ifelse(positive.Al2==1, "qPCR +", "qPCR -"),
         positive.Hw.lab = ifelse(positive.Hw==1, "qPCR +", "qPCR -"),
         positive.Tt.lab = ifelse(positive.Tt==1, "qPCR +", "qPCR -"),
         alkk.lab = ifelse(alkk==1, "KK +", "KK -"),
         hwkk.lab = ifelse(hwkk==1, "KK +", "KK -"),
         ttkk.lab = ifelse(ttkk==1, "KK +", "KK -"))

qdata = qdata %>% 
  mutate(positive.Al2.lab = factor(positive.Al2.lab, levels=c("qPCR +", "qPCR -")),
         positive.Hw.lab = factor(positive.Hw.lab, levels=c("qPCR +", "qPCR -")),
         positive.Tt.lab = factor(positive.Tt.lab, levels=c("qPCR +", "qPCR -")),
         alkk.lab = factor(alkk.lab, levels=c("KK +", "KK -")),
         hwkk.lab = factor(hwkk.lab, levels=c("KK +", "KK -")),
         ttkk.lab = factor(ttkk.lab, levels=c("KK +", "KK -")))

#---------------------------------------
# create table
#---------------------------------------
altab = table(qdata$alkk.lab, qdata$positive.Al2.lab)
hwtab = table(qdata$hwkk.lab, qdata$positive.Hw.lab)
tttab = table(qdata$ttkk.lab, qdata$positive.Tt.lab)

alptab = prop.table(table(qdata$alkk.lab, qdata$positive.Al2.lab))
hwptab = prop.table(table(qdata$hwkk.lab, qdata$positive.Hw.lab))
ttptab = prop.table(table(qdata$ttkk.lab, qdata$positive.Tt.lab))

al_results = make2x2(altab, alptab, label = "A. lumbricoides")
hw_results = make2x2(hwtab, hwptab, label = "Hookworm")
tt_results = make2x2(tttab, ttptab, label = "T. trichiura")

all_tables = bind_rows(al_results, hw_results, tt_results)

all_tables = all_tables %>% 
  mutate(kappa = c("", paste0(round(al.kappa$value, 2), " (<0.001)"), "",
                   "", paste0(round(hw.kappa$value, 2), " (<0.001)"), "",
                   "", paste0(round(tt.kappa$value, 2), " (<0.001)"), ""))

#---------------------------------------
# save table
#---------------------------------------
write.csv(all_tables, file=paste0(tab_dir, "2x2.table.csv"),row.names=FALSE)
