#######################################
# WASH Benefits Bangladesh STH KK qPCR validation
# estimate sensitivity, specificity
# of kk using qPCR as gold standard
# make table
#######################################
rm(list=ls())
library(dplyr)

load("~/Dropbox/WASH-B-STH-Add-on/TFGH/Data/RData/sensspec_qgold.RData")

source("~/Documents/CRG/wash-benefits/bangladesh/src/wbb-sth-qpcr/2-fig-tab/0-base-table-functions.R")


sens=kk.sens.qgold %>%
  rename(lb=Lower.95.CI, ub=Upper.95.CI) %>%
  mutate(sens=ptestci.format(Mean,lb,ub,decimals=1,scale=100)) %>%
  select(org,N,sens)

spec=kk.spec.qgold %>%
  rename(lb=Lower.95.CI, ub=Upper.95.CI) %>%
  mutate(spec=ptestci.format(Mean,lb,ub,decimals=1,scale=100)) %>%
  select(org,N,spec)

out=full_join(sens,spec,by=c("org","N"))

write.csv(out,file="~/Dropbox/WASH-B-STH-Add-on/TFGH/Results/ss_qgold.csv",row.names=FALSE)
