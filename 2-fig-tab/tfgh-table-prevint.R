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

source("~/Documents/CRG/wash-benefits/bangladesh/src/wbb-sth-qpcr/2-fig-tab/0-base-table-functions.R")

prev.kk=as.data.frame(rbind(al.kk,hw.kk,tt.kk))
prev.kk$org=c("Ascaris lumbricoides","Hookworm","Trichuris trichiura")
colnames(prev.kk)[5:6]=c("lb","ub")
prev.kk = prev.kk %>%
  mutate(prevkk=ptestci.format(Mean,lb,ub,decimals=1,scale=100)) %>%
  select(org,N,prevkk) %>%
  rename(Nkk=N)

prev.q=as.data.frame(rbind(al.q,na.q,ac.q,ad.q,tt.q))
prev.q$org=c("Ascaris lumbricoides","Trichuris trichiura",
           "Necator americanus","Ancylostoma ceylanicum","Ancylostoma duodenale")
colnames(prev.q)[5:6]=c("lb","ub")
prev.q = prev.q %>%
  mutate(prevq=ptestci.format(Mean,lb,ub,decimals=1,scale=100)) %>%
  select(org,N,prevq)%>%
  rename(Nq=N)

epg=as.data.frame(rbind(al.kk.gmn,hw.kk.gmn,tt.kk.gmn))
epg$org=c("Ascaris lumbricoides","Hookworm","Trichuris trichiura")
colnames(epg)[5:6]=c("lb","ub")
epg = epg %>%
  mutate(epg=ptestci.format(Mean,lb,ub,decimals=1,scale=100)) %>%
  select(org,epg)

dna=as.data.frame(rbind(al.q.gmn,na.q.gmn,ac.q.gmn,ad.q.gmn,tt.q.gmn))
dna$org=c("Ascaris lumbricoides","Trichuris trichiura",
          "Necator americanus","Ancylostoma ceylanicum","Ancylostoma duodenale")
colnames(dna)[5:6]=c("lb","ub")
dna = dna %>%
  mutate(dna=ptestci.format(Mean,lb,ub,decimals=1,scale=100)) %>%
  select(org,dna)

tab=full_join(prev.kk,prev.q,by="org") 
tab=full_join(tab,epg,by="org")
tab=full_join(tab,dna,by="org")

tab = tab %>% 
  select(org,Nkk,prevkk,epg,Nq,prevq,dna) %>%
  mutate_at(.vars=vars(Nkk,prevkk,Nq,prevq,epg,dna),
            .funs=repNA) 

# manual reorder
tab=tab[c(1,2,4:6,3),]

write.csv(tab,file="~/Dropbox/WASH-B-STH-Add-on/TFGH/Results/prev.table.csv",row.names=FALSE)
