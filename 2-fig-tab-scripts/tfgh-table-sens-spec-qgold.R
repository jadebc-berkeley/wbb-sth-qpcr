#######################################
# WASH Benefits Bangladesh STH KK qPCR validation
# estimate sensitivity, specificity
# of kk using qPCR as gold standard
# make table

# by Jade Benjamin-Chung
# jadebc@berkeley.edu
#######################################
rm(list=ls())

# configure directories, load libraries and base functions
source(paste0(here::here(), "/0-config.R"))

# load data
load(paste0(data_dir,"sensspec_qgold.RData"))
load(paste0(data_dir,"qdata.RData"))

#---------------------------------------
# make table
#---------------------------------------
N.al=nrow(qdata[!is.na(qdata$positive.Al) & !is.na(qdata$alkk),])
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

#---------------------------------------
# save table
#---------------------------------------
write.csv(out,file=paste0(tab_dir, "ss_qgold.csv"),row.names=FALSE)
