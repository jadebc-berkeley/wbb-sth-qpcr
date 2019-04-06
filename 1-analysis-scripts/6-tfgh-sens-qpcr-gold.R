#######################################
# WASH Benefits Bangladesh STH KK qPCR validation
# estimate sensitivity, specificity
# of kk using qPCR as gold standard

# by Jade Benjamin-Chung
# jadebc@berkeley.edu
#######################################
rm(list=ls())

# configure directories, load libraries and base functions
source(paste0(here::here(), "/0-config.R"))

# load data
load(paste0(data_dir,"qdata.RData"))

#----------------------------------------
# subset to relevant columns
#----------------------------------------
d = qdata %>%
  select(dataid,personid,clusterid,positive.Hw,
         positive.Tt,alkk,hwkk,ttkk)

#----------------------------------------
# calculate sensitivity
#----------------------------------------
ald1= d %>% filter(positive.Al==1) %>%
  filter(!is.na(alkk)) 
al.sens=washb_mean(ald1$alkk,id=ald1$clusterid)

hwd1= d %>% filter(positive.Hw==1) %>%
  filter(!is.na(hwkk)) 
hw.sens=washb_mean(hwd1$hwkk,id=hwd1$clusterid)

ttd1= d %>% filter(positive.Tt==1) %>%
  filter(!is.na(ttkk)) 
tt.sens=washb_mean(ttd1$alkk,id=ttd1$clusterid)

kk.sens.qgold=data.frame(rbind(al.sens,hw.sens,tt.sens))
kk.sens.qgold$org=c("A. lumbricoides","Hookworm","T. trichiura")

#----------------------------------------
# calculate specificity
#----------------------------------------
ald2= d %>% filter(positive.Al==0) %>%
  filter(!is.na(alkk)) %>%
  mutate(alneg=ifelse(alkk==1,0,1))
al.spec=washb_mean(ald2$alneg,id=ald2$clusterid)

hwd2= d %>% filter(positive.Hw==0) %>%
  filter(!is.na(hwkk)) %>%
  mutate(hwneg=ifelse(hwkk==1,0,1))
hw.spec=washb_mean(hwd2$hwneg,id=hwd2$clusterid)

ttd2= d %>% filter(positive.Tt==0) %>%
  filter(!is.na(ttkk)) %>%
  mutate(ttneg=ifelse(ttkk==1,0,1))
tt.spec=washb_mean(ttd2$ttneg,id=ttd2$clusterid)

kk.spec.qgold=data.frame(rbind(al.spec,hw.spec,tt.spec))
kk.spec.qgold$org=c("A. lumbricoides","Hookworm","T. trichiura")

save(kk.sens.qgold, kk.spec.qgold,
     file=paste0(data_dir, "sensspec_qgold.RData"))

