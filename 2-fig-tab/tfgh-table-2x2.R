#######################################
# WASH Benefits Bangladesh STH KK qPCR validation
# table of concordance discordance
#######################################
rm(list=ls())
library(dplyr)
library(tidyr)
library(reshape2)
library(ggplot2)
library(grid)
library(gridExtra)
load("~/Dropbox/WASH-B-STH-Add-on/TFGH/Data/RData/qdata.RData")

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

altab = table(qdata$alkk.lab, qdata$positive.Al2.lab)
hwtab = table(qdata$hwkk.lab, qdata$positive.Hw.lab)
tttab = table(qdata$ttkk.lab, qdata$positive.Tt.lab)

alptab = prop.table(table(qdata$alkk.lab, qdata$positive.Al2.lab))
hwptab = prop.table(table(qdata$hwkk.lab, qdata$positive.Hw.lab))
ttptab = prop.table(table(qdata$ttkk.lab, qdata$positive.Tt.lab))