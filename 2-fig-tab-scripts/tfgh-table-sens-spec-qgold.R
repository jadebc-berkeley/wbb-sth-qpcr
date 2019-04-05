#######################################
# WASH Benefits Bangladesh STH KK qPCR validation
# estimate sensitivity, specificity
# of kk using qPCR as gold standard
# make table
#######################################
rm(list=ls())
library(dplyr)

load("~/Dropbox/WASH-B-STH-Add-on/TFGH/Data/RData/sensspec_qgold.RData")
load("~/Dropbox/WASH-B-STH-Add-on/TFGH/Data/RData/qdata.RData")

source("~/Documents/CRG/wash-benefits/bangladesh/src/wbb-sth-qpcr/2-fig-tab/0-base-table-functions.R")

N.al=nrow(qdata[!is.na(qdata$positive.Al2) & !is.na(qdata$alkk),])
N.hw=nrow(qdata[!is.na(qdata$positive.Hw) & !is.na(qdata$hwkk),])
N.tt=nrow(qdata[!is.na(qdata$positive.Tt) & !is.na(qdata$ttkk),])

sens=kk.sens.qgold %>%
  select(org, Mean, lb=Lower.95.CI, ub=Upper.95.CI) %>%
  mutate(sens=ptestci.format(Mean,lb,ub,decimals=0,scale=100)) %>%
  select(org,sens)

spec=kk.spec.qgold %>%
  select(org, Mean, lb=Lower.95.CI, ub=Upper.95.CI) %>%
  mutate(spec=ptestci.format(Mean,lb,ub,decimals=0,scale=100)) %>%
  select(org,spec)

out=full_join(sens,spec,by=c("org"))
out$N=c(N.al,N.hw,N.tt)

write.csv(out,file="~/Dropbox/WASH-B-STH-Add-on/TFGH/Results/ss_qgold.csv",row.names=FALSE)
