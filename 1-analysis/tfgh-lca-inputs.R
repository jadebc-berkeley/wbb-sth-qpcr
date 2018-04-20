#######################################
# WASH Benefits Bangladesh STH KK qPCR validation
# inputs for bayesian lca
#######################################
rm(list=ls())
library(dplyr)
load("~/Dropbox/WASH-B-STH-Add-on/TFGH/Data/RData/qdata.RData")

#--------------------------------------
# hookworm
#--------------------------------------
qdata=qdata[!is.na(qdata$positive.Hw) & !is.na(qdata$hwkk),]
nrow(qdata)

nrow(qdata[qdata$positive.Hw==1 & qdata$hwkk==1,])
nrow(qdata[qdata$positive.Hw==0 & qdata$hwkk==1,])
nrow(qdata[qdata$positive.Hw==1 & qdata$hwkk==0,])
nrow(qdata[qdata$positive.Hw==0 & qdata$hwkk==0,])

# initial value for prevalence - use qpcr
mean(qdata$positive.Hw)


