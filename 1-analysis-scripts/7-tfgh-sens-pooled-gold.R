#######################################
# WASH Benefits Bangladesh STH KK qPCR validation
# estimate sensitivity, specificity
# of kk using either method as gold

# by Jade Benjamin-Chung
# jadebc@berkeley.edu
#######################################
rm(list=ls())

# configure directories, load libraries and base functions
source(paste0(here::here(), "/0-config.R"))

# load data
load(paste0(data_dir,"qdata.RData"))

#--------------------------------------
# subset to relevant columns
#--------------------------------------
d = qdata %>%
  select(dataid,personid,clusterid,positive.Al,positive.Hw,
         positive.Tt,alkk,hwkk,ttkk)

#----------------------------------------
# calculate sensitivity
#----------------------------------------
ald1= d %>% filter(positive.Al==1) %>%
  filter(!is.na(alkk)) %>%
  filter(!is.na(positive.Al)) 
al.kk.sens=washb_mean(ald1$alkk,id=ald1$clusterid)
al.q.sens=washb_mean(ald1$positive.Al,id=ald1$clusterid)

hwd1= d %>% filter(hwkk==1 | positive.Hw==1) %>%
  filter(!is.na(hwkk)) %>%
  filter(!is.na(positive.Hw)) 
hw.kk.sens=washb_mean(hwd1$hwkk,id=hwd1$clusterid)
hw.q.sens=washb_mean(hwd1$positive.Hw,id=hwd1$clusterid)

ttd1= d %>% filter(ttkk==1 | positive.Tt==1) %>%
  filter(!is.na(ttkk)) %>%
  filter(!is.na(positive.Tt)) 
tt.kk.sens=washb_mean(ttd1$ttkk,id=ttd1$clusterid)
tt.q.sens=washb_mean(ttd1$positive.Tt,id=ttd1$clusterid)

kk.sens.pooled.gold=data.frame(rbind(al.kk.sens,hw.kk.sens,tt.kk.sens))
kk.sens.pooled.gold$org=c("A. lumbricoides","Hookworm","T. trichiura")

q.sens.pooled.gold=data.frame(rbind(al.q.sens,hw.q.sens,tt.q.sens))
q.sens.pooled.gold$org=c("A. lumbricoides","Hookworm","T. trichiura")

#----------------------------------------
# calculate specificity
#----------------------------------------
ald2= d %>% filter(positive.Al==0) %>%
  filter(!is.na(alkk)) %>%
  filter(!is.na(positive.Al)) %>% 
  mutate(alkkneg=ifelse(alkk==1,0,1),
         alqneg = ifelse(positive.Al==1, 0, 1))
al.kk.spec=washb_mean(ald2$alkkneg,id=ald2$clusterid)
al.q.spec=washb_mean(ald2$alqneg,id=ald2$clusterid)

hwd2= d %>% filter(positive.Hw==0 | hwkk==0) %>%
  filter(!is.na(hwkk)) %>%
  filter(!is.na(positive.Hw)) %>%
  mutate(hwkkneg=ifelse(hwkk==1,0,1),
         hwqneg = ifelse(positive.Hw==1,0,1))
hw.kk.spec=washb_mean(hwd2$hwkkneg,id=hwd2$clusterid)
hw.q.spec=washb_mean(hwd2$hwqneg,id=hwd2$clusterid)

ttd2= d %>% filter(positive.Tt==0 | ttkk==0) %>%
  filter(!is.na(ttkk)) %>%
  filter(!is.na(positive.Tt)) %>%
  mutate(ttkkneg=ifelse(ttkk==1,0,1),
         ttqneg=ifelse(positive.Tt==1,0,1))
tt.kk.spec=washb_mean(ttd2$ttkkneg,id=ttd2$clusterid)
tt.q.spec=washb_mean(ttd2$ttqneg,id=ttd2$clusterid)

kk.spec.pooled.gold=data.frame(rbind(al.kk.spec,hw.kk.spec,tt.kk.spec))
kk.spec.pooled.gold$org=c("A. lumbricoides","Hookworm","T. trichiura")

q.spec.pooled.gold=data.frame(rbind(al.q.spec,hw.q.spec,tt.q.spec))
q.spec.pooled.gold$org=c("A. lumbricoides","Hookworm","T. trichiura")

#--------------------------------------
# save results
#--------------------------------------
save(kk.sens.pooled.gold, kk.spec.pooled.gold,
     q.sens.pooled.gold, q.spec.pooled.gold,
     file=paste0(data_dir,"sensspec_pooled_gold.RData"))

