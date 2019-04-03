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


#--------------------------------------
# ascaris
#--------------------------------------
qdata=qdata[!is.na(qdata$positive.Al) & !is.na(qdata$alkk),]
nrow(qdata)

nrow(qdata[qdata$positive.Al==1 & qdata$alkk==1,])
nrow(qdata[qdata$positive.Al==0 & qdata$alkk==1,])
nrow(qdata[qdata$positive.Al==1 & qdata$alkk==0,])
nrow(qdata[qdata$positive.Al==0 & qdata$alkk==0,])

# initial value for prevalence - use qpcr
mean(qdata$positive.Al)

# new assay
qdata=qdata[!is.na(qdata$positive.Al2) & !is.na(qdata$alkk),]
nrow(qdata)

nrow(qdata[qdata$positive.Al2==1 & qdata$alkk==1,])
nrow(qdata[qdata$positive.Al2==0 & qdata$alkk==1,])
nrow(qdata[qdata$positive.Al2==1 & qdata$alkk==0,])
nrow(qdata[qdata$positive.Al2==0 & qdata$alkk==0,])

# initial value for prevalence - use qpcr
mean(qdata$positive.Al2)


#--------------------------------------
# trichuris
#--------------------------------------
qdata=qdata[!is.na(qdata$positive.Tt) & !is.na(qdata$ttkk),]
nrow(qdata)

nrow(qdata[qdata$positive.Tt==1 & qdata$ttkk==1,])
nrow(qdata[qdata$positive.Tt==0 & qdata$ttkk==1,])
nrow(qdata[qdata$positive.Tt==1 & qdata$ttkk==0,])
nrow(qdata[qdata$positive.Tt==0 & qdata$ttkk==0,])

# initial value for prevalence - use qpcr
mean(qdata$positive.Tt)




