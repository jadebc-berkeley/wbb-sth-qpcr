rm(list=ls())
load("~/Dropbox/WASH-B-STH-Add-on/TFGH/Data/RData/qdata.RData")

prop.table(table(qdata$sex))

qdata %>%
  summarise(min=min(agem,na.rm=TRUE),
            max=max(agey,na.rm=TRUE),
            mean=mean(agem,na.rm=TRUE))
