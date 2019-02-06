#######################################
# WASH Benefits Bangladesh STH KK qPCR validation
# characteristics of study participants
#######################################
rm(list=ls())
library(dplyr)

load("~/Dropbox/WASH-B-STH-Add-on/TFGH/Data/RData/qdata.RData")

prop.table(table(qdata$sex))
prop.table(table(qdata$dw))

qdata %>%
  summarise(min=min(agem,na.rm=TRUE),
            max=max(agem,na.rm=TRUE),
            mean=mean(agem,na.rm=TRUE))
