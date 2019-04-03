#######################################
# WASH Benefits Bangladesh STH KK qPCR validation
# Subset ids requested by Smith for
# further analysis
#######################################
rm(list=ls())
library(reshape2)
library(ggplot2)
library(dplyr)
library(tidyr)

load("~/Dropbox/WASH-B-STH-Add-on/TFGH/Data/RData/qdata.RData")

#--------------------------------------
# ids with positive ascaris kk and 
# all negative for qPCR
#--------------------------------------
qdata$keep=ifelse(qdata$alepg>0 & qdata$positive.Ad==0 &
                    qdata$positive.Na==0 & qdata$positive.Tt==0 & 
                    qdata$positive.Al==0 & qdata$positive.Ac==0 &
                    qdata$hwepg==0 & qdata$ttepg==0,1,0)

keep=qdata[qdata$keep==1,]
keep=keep[!is.na(keep$alint),]
keep=keep[rev(order(keep$alepg)),]
keep[1:6,c("dataid","personid")]
keep[7:16,c("dataid","personid")]
keep[17:29,c("dataid","personid")]

# requested by Nils on 2/22/18
# kk and qpcr positive for tt and nothing else
qdata$keep2=ifelse(qdata$ttkk==1 & qdata$hwkk==0 & qdata$alkk==0 & 
                     qdata$positive.Tt==1 & qdata$positive.Na!=1 & qdata$positive.Ad!=1 & 
                     qdata$positive.Ss!=1 & qdata$positive.Al!=1 & qdata$positive.Ac!=1,1,0)
as.matrix(paste0(qdata$dataid[qdata$keep2==1],"E",qdata$personid[qdata$keep2==1],"S1"))

# kk and qpcr positive for Na and nothing else
qdata$keep3=ifelse(qdata$ttkk==0 & qdata$hwkk==1 & qdata$alkk==0 & 
                     qdata$positive.Tt!=1 & qdata$positive.Na==1 & qdata$positive.Ad!=1 & 
                     qdata$positive.Ss!=1 & qdata$positive.Al!=1 & qdata$positive.Ac!=1,1,0)
as.matrix(paste0(qdata$dataid[qdata$keep3==1],"E",qdata$personid[qdata$keep3==1],"S1"))

# kk and qpcr positive for Al and nothing else
qdata$keep4=ifelse(qdata$ttkk==0 & qdata$hwkk==0 & qdata$alkk==1 & 
                     qdata$positive.Tt!=1 & qdata$positive.Na!=1 & qdata$positive.Ad!=1 & 
                     qdata$positive.Ss!=1 & qdata$positive.Al==1 & qdata$positive.Ac!=1,1,0)
as.matrix(paste0(qdata$dataid[qdata$keep4==1],"E",qdata$personid[qdata$keep4==1],"S1"))

# kk and qpcr negative for all
qdata$keep5=ifelse(qdata$ttkk==0 & qdata$hwkk==0 & qdata$alkk==0 & 
                     qdata$positive.Tt!=1 & qdata$positive.Na!=1 & qdata$positive.Ad!=1 & 
                     qdata$positive.Ss!=1 & qdata$positive.Al!=1 & qdata$positive.Ac!=1,1,0)
as.matrix(paste0(qdata$dataid[qdata$keep5==1],"E",qdata$personid[qdata$keep5==1],"S1")[1:15])


