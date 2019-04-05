#######################################
# WASH Benefits Bangladesh STH KK qPCR validation
# table of estimates of sensitivity, specificity
# using pooled kk qPCR as gold standard
#######################################
rm(list=ls())
library(dplyr)

load("~/Dropbox/WASH-B-STH-Add-on/TFGH/Data/RData/sensspec_pooled_gold.RData")
load("~/Dropbox/WASH-B-STH-Add-on/TFGH/Data/RData/qdata.RData")

source("~/Documents/CRG/wash-benefits/bangladesh/src/wbb-sth-qpcr/2-fig-tab/0-base-table-functions.R")

N.al=nrow(qdata[!is.na(qdata$positive.Al2) & !is.na(qdata$alkk),])
N.hw=nrow(qdata[!is.na(qdata$positive.Hw) & !is.na(qdata$hwkk),])
N.tt=nrow(qdata[!is.na(qdata$positive.Tt) & !is.na(qdata$ttkk),])

kk.sens=kk.sens.pooled.gold %>%
  select(org, Mean, lb=Lower.95.CI, ub=Upper.95.CI) %>%
  mutate(sens=ptestci.format(Mean,lb,ub,decimals=0,scale=100)) %>%
  select(org,sens)

kk.spec=kk.spec.pooled.gold %>%
  select(org, Mean, lb=Lower.95.CI, ub=Upper.95.CI) %>%
  mutate(spec=ptestci.format(Mean,lb,ub,decimals=0,scale=100)) %>%
  select(org,spec)

q.sens=q.sens.pooled.gold %>%
  select(org, Mean, lb=Lower.95.CI, ub=Upper.95.CI) %>%
  mutate(sens=ptestci.format(Mean,lb,ub,decimals=0,scale=100)) %>%
  select(org,sens)

q.spec=q.spec.pooled.gold %>%
  select(org, Mean, lb=Lower.95.CI, ub=Upper.95.CI) %>%
  mutate(spec=ptestci.format(Mean,lb,ub,decimals=0,scale=100)) %>%
  select(org,spec)

out=full_join(kk.sens,kk.spec,by=c("org"))
out=full_join(out,q.sens,by=c("org"))
out=full_join(out,q.spec,by=c("org"))
out$N=c(N.al,N.hw,N.tt)

write.csv(out,file="~/Dropbox/WASH-B-STH-Add-on/TFGH/Results/ss_pooled_gold.csv",row.names=FALSE)
