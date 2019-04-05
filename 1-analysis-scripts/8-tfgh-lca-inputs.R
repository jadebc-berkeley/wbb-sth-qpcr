#######################################
# WASH Benefits Bangladesh STH KK qPCR validation
# produce initial values for bayesian lca models
#######################################
rm(list=ls())

# configure directories, load libraries and base functions
source(paste0(here::here(), "/0-config.R"))

# load data
load(paste0(data_dir,"qdata.RData"))

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




