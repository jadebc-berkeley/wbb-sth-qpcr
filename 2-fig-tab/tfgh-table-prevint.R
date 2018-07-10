#######################################
# WASH Benefits Bangladesh STH KK qPCR validation
# table of qpcr vs. kk prevalence and intensity
#######################################
rm(list=ls())
library(dplyr)
library(tidyr)
library(reshape2)
library(ggplot2)
library(gridExtra)
load("~/Dropbox/WASH-B-STH-Add-on/TFGH/Data/RData/qdata.RData")
load("~/Dropbox/WASH-B-STH-Add-on/TFGH/Data/RData/prev_results.RData")
load("~/Dropbox/WASH-B-STH-Add-on/TFGH/Data/RData/concentration.RData")

source("~/Documents/CRG/wash-benefits/bangladesh/src/wbb-sth-qpcr/2-fig-tab/0-base-table-functions.R")

# -------------------------------------
# prevalence
# -------------------------------------
prev.kk=as.data.frame(rbind(al.kk,hw.kk,tt.kk))
prev.kk$org=c("Ascaris lumbricoides","Hookworm","Trichuris trichiura")
colnames(prev.kk)[5:6]=c("lb","ub")
prev.kk = prev.kk %>%
  mutate(prevkk=ptestci.format(Mean,lb,ub,decimals=1,scale=100)) %>%
  select(org,Nkk=N,prevkk) 

prev.q=as.data.frame(rbind(al.q,hw.q,na.q,ac.q,ad.q,tt.q))
prev.q$org=c("Ascaris lumbricoides","Hookworm",
           "Necator americanus","Ancylostoma ceylanicum","Ancylostoma duodenale",
           "Trichuris trichiura")
colnames(prev.q)[5:6]=c("lb","ub")
prev.q = prev.q %>%
  mutate(prevq=ptestci.format(Mean,lb,ub,decimals=1,scale=100)) %>%
  select(org,Nq=N,prevq)
# -------------------------------------
# KK EPG
# -------------------------------------
epg=as.data.frame(rbind(al.kk.gmn,hw.kk.gmn,tt.kk.gmn))
epg$org=c("Ascaris lumbricoides","Hookworm","Trichuris trichiura")
colnames(epg)[5:6]=c("lb","ub")
epg = epg %>%
  mutate(epg=ptestci.format(Mean,lb,ub,decimals=1,scale=100)) %>%
  select(org,epg)

# -------------------------------------
# DNA concentration
# -------------------------------------
# function to make pretty median (range)
medrange=function(y){
  min=quantile(y,prob=c(0))
  med=quantile(y,prob=c(0.5))
  max=quantile(y,prob=c(1))
  
  min=formatC(min, format = "e", digits = 2)
  # med=formatC(med, format = "e", digits = 2)
  med=sprintf("%0.1f",med)
  # max=sprintf("%0.1f",as.numeric(max))
  max=formatC(max, format = "e", digits = 2)
  return(paste0(med, " (",min, ", ",max,")"))
}

al.dna=medrange(al$copies)
na.dna=medrange(na$copies)
ac.dna=medrange(ac$copies)
ad.dna=medrange(ad$copies)
tt.dna=medrange(tt$copies)

dna=as.data.frame(rbind(al.dna,na.dna,ac.dna,ad.dna,tt.dna))
org=c("Ascaris lumbricoides","Trichuris trichiura",
          "Necator americanus","Ancylostoma ceylanicum","Ancylostoma duodenale")
dna=cbind(org,dna)
dna$org=as.character(dna$org)
colnames(dna)[2]="conc"

tab=full_join(prev.kk,prev.q,by="org") 
tab=full_join(tab,epg,by="org")
tab=full_join(tab,dna,by="org")

# manual reorder
tab = tab[,c(1:3,6,4,5,7)]
tab=tab[c(1,2,4:6,3),]

write.csv(tab,file="~/Dropbox/WASH-B-STH-Add-on/TFGH/Results/prev.table.csv",row.names=FALSE)
