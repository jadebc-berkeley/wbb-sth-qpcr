#######################################
# WASH Benefits Bangladesh STH KK qPCR validation
# characteristics of study participants
#######################################
rm(list=ls())

# configure directories, load libraries and base functions
source(paste0(here::here(), "/0-config.R"))

# load data
load(paste0(data_dir,"qdata.RData"))


prop.table(table(qdata$sex))
prop.table(table(qdata$dw))

qdata %>%
  summarise(min=min(agem,na.rm=TRUE),
            max=max(agem,na.rm=TRUE),
            mean=mean(agem,na.rm=TRUE))
