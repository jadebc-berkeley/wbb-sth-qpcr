#######################################
# WASH Benefits Bangladesh STH KK qPCR validation
# calculate prevalence and geometric mean
# using each diagnostic
# 
#######################################
rm(list=ls())
library(dplyr)
library(clusrank)
library(tidyr)
library(washb)
load("~/Dropbox/WASH-B-STH-Add-on/TFGH/Data/RData/qdata.RData")
load("~/Dropbox/WASH-B-STH-Add-on/TFGH/Data/RData/concentration.RData")

#--------------------------------------
# Estimate prevalence and 95%CI for each
# test, robust sandwich SEs
#--------------------------------------
al.kk=washb_mean(qdata$alkk, id=qdata$clusterid, print = TRUE)
hw.kk=washb_mean(qdata$hwkk, id=qdata$clusterid, print = TRUE)
tt.kk=washb_mean(qdata$ttkk, id=qdata$clusterid, print = TRUE)

al.q=washb_mean(qdata$positive.Al, id=qdata$clusterid, print = TRUE)
hw.q=washb_mean(qdata$positive.Hw, id=qdata$clusterid, print = TRUE)
na.q=washb_mean(qdata$positive.Na, id=qdata$clusterid, print = TRUE)
ac.q=washb_mean(qdata$positive.Ac, id=qdata$clusterid, print = TRUE)
ad.q=washb_mean(qdata$positive.Ad, id=qdata$clusterid, print = TRUE)
tt.q=washb_mean(qdata$positive.Tt, id=qdata$clusterid, print = TRUE)
ss.q=washb_mean(qdata$positive.Ss, id=qdata$clusterid, print = TRUE)

#--------------------------------------
# co-infection
#--------------------------------------
qdata <- qdata %>%
  mutate(numsth.kk=alkk+hwkk+ttkk,
         numsth.q=positive.Al+positive.Ac+positive.Na+
           positive.Ad+positive.Tt) %>%
  mutate(multisth.kk=ifelse(numsth.kk>1,1,0),
         multisth.q=ifelse(numsth.q>1,1,0))

prop.table(table(qdata$numsth.kk))*100
prop.table(table(qdata$numsth.q))*100

prop.table(table(qdata$multisth.kk))*100
prop.table(table(qdata$multisth.q))*100

#--------------------------------------
# mod/heavy intensity infection
#--------------------------------------
qdata.conc <- qdata.conc %>%
  mutate(almh=ifelse(alepg>=5000,1,0),
         hwmh=ifelse(hwepg>=2000,1,0),
         ttmh=ifelse(ttepg>=1000,1,0),
         almh.f=as.factor(ifelse(almh==1,"Moderate-heavy intensity\ninfection","Low intensity\ninfection")),
         hwmh.f=as.factor(ifelse(hwmh==1,"Moderate-heavy intensity","Low intensity")),
         ttmh.f=as.factor(ifelse(ttmh==1,"Moderate-heavy intensity","Low intensity")))

prop.table(table(qdata.conc$almh))*100
prop.table(table(qdata.conc$hwmh))*100
prop.table(table(qdata.conc$ttmh))*100

quantile(qdata.conc$CTmean.Al[qdata.conc$almh==1],probs=c(0,0.5,1),na.rm=TRUE)
quantile(qdata.conc$CTmean.Al[qdata.conc$almh==0],probs=c(0,0.5,1),na.rm=TRUE)

pdf(file="~/Dropbox/WASH-B-STH-Add-on/TFGH/Results/wbb-alconc-boxplot.pdf",
    width=5,height=5)
ggplot(qdata.conc[!is.na(qdata.conc$almh.f),],aes(y=log10(copies.Al),x=almh.f))+
  geom_boxplot()+  
  geom_dotplot(aes(fill=almh.f,col=almh.f),binaxis='y',stackdir='center',
               stackratio=1.3,dotsize=0.7,binwidth=.15, alpha=0.5)+
  xlab("Infection intensity")+
  scale_color_manual(values=c("#CB59EB","#E37F2D"),guide=FALSE)+
  scale_fill_manual(values=c("#CB59EB","#E37F2D"),guide=FALSE)+
  ylab(expression(paste("Mean", " log"[10], italic(" A. lumbricoides"), " DNA (ag/",mu,"l)")))+
  theme_bw()
dev.off()

qdata.conc <- qdata.conc %>%
  mutate(almh.alt=ifelse(alepg>=15000,1,0),
         almh.alt.f=as.factor(ifelse(almh.alt==1,"Moderate-heavy intensity\ninfection","Low intensity\ninfection")))

pdf(file="~/Dropbox/WASH-B-STH-Add-on/TFGH/Results/wbb-alconc-boxplot.pdf",
    width=5,height=5)
ggplot(qdata.conc[!is.na(qdata.conc$almh.alt.f),],aes(y=log10(copies.Al),x=almh.alt.f))+
  geom_boxplot()+  
  geom_dotplot(aes(fill=almh.alt.f,col=almh.alt.f),binaxis='y',stackdir='center',
               stackratio=1.3,dotsize=0.7,binwidth=.15, alpha=0.5)+
  xlab("Infection intensity")+
  scale_color_manual(values=c("#CB59EB","#E37F2D"),guide=FALSE)+
  scale_fill_manual(values=c("#CB59EB","#E37F2D"),guide=FALSE)+
  ylab(expression(paste("Mean", " log"[10], italic(" A. lumbricoides"), " DNA (ag/",mu,"l)")))+
  theme_bw()
dev.off()

#--------------------------------------
# Estimate geometric mean and 95%CI for each
# test, robust sandwich SEs
#--------------------------------------
qdata = qdata %>%
  mutate(ln.alepg=log(alepg+1),
         ln.hwepg=log(hwepg+1),
         ln.ttepg=log(ttepg+1),
         ln.CT.al=log(CTmean.Al+1),
         ln.CT.na=log(CTmean.Na+1),
         ln.CT.ac=log(CTmean.Ac+1),
         ln.CT.ad=log(CTmean.Ad+1),
         ln.CT.tt=log(CTmean.Tt+1))

al.kk.gmn=washb_mean(qdata$ln.alepg, id=qdata$clusterid, print = TRUE)
hw.kk.gmn=washb_mean(qdata$ln.hwepg, id=qdata$clusterid, print = TRUE)
tt.kk.gmn=washb_mean(qdata$ln.ttepg, id=qdata$clusterid, print = TRUE)

al.q.gmn=washb_mean(qdata$ln.CT.al, id=qdata$clusterid, print = TRUE)
na.q.gmn=washb_mean(qdata$ln.CT.na, id=qdata$clusterid, print = TRUE)
ac.q.gmn=washb_mean(qdata$ln.CT.ac, id=qdata$clusterid, print = TRUE)
ad.q.gmn=washb_mean(qdata$ln.CT.ad, id=qdata$clusterid, print = TRUE)
tt.q.gmn=washb_mean(qdata$ln.CT.tt, id=qdata$clusterid, print = TRUE)


save(al.kk,hw.kk,tt.kk,
     al.q,hw.q,na.q,ac.q,ad.q,tt.q,
     al.kk.gmn,hw.kk.gmn,tt.kk.gmn,
     al.q.gmn,na.q.gmn,ac.q.gmn,ad.q.gmn,tt.q.gmn,
     file="~/Dropbox/WASH-B-STH-Add-on/TFGH/Data/RData/prev_results.RData")
