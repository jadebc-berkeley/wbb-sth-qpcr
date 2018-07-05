#######################################
# WASH Benefits Bangladesh STH KK qPCR validation
# sensitivity and specificity analysis
#######################################
rm(list=ls())
library(dplyr)
library(clusrank)
library(tidyr)
library(washb)
load("~/Dropbox/WASH-B-STH-Add-on/TFGH/Data/RData/qdata.RData")


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


#--------------------------------------
# Compare prevalence using 
# Wilcoxon matched-pairs signed rank test 
# appropriate for non-normally distributed outcomes

# assumes data are independent.. but they are not
#--------------------------------------
qdata$personidn=ifelse(qdata$personid=="C1",1,
                       ifelse(qdata$personid=="E1",2,
                              ifelse(qdata$personid=="O1",3,
                                     ifelse(qdata$personid=="T1",4,5))))
# prepare data
al.dat=qdata %>%
  select(dataid,clusterid,personid,personidn,alkk,positive.Al) %>%
  gather(test,positive, 5:6) %>%
  mutate(test=ifelse(test=="alkk","Kato-Katz","qPCR")) %>%
  mutate(id=as.numeric(paste0(dataid,personidn))) %>%
  arrange(id,test) %>% select(id,positive,test)

hw.dat=qdata %>%
  select(dataid,clusterid,personid,personidn,hwkk,positive.Hw) %>%
  gather(test,positive, 5:6) %>%
  mutate(test=ifelse(test=="hwkk","Kato-Katz","qPCR")) %>%
  mutate(id=as.numeric(paste0(dataid,personidn))) %>%
  arrange(id,test) %>% select(id,positive,test)

tt.dat=qdata %>%
  select(dataid,clusterid,personid,personidn,ttkk,positive.Tt) %>%
  gather(test,positive, 5:6) %>%
  mutate(test=ifelse(test=="ttkk","Kato-Katz","qPCR")) %>%
  mutate(id=as.numeric(paste0(dataid,personidn))) %>%
  arrange(id,test) %>% select(id,positive,test)

al.wilcox=clusWilcox.test(positive ~ test,data=al.dat,paired=TRUE)
hw.wilcox=clusWilcox.test(positive ~ test,data=hw.dat,paired=TRUE)
tt.wilcox=clusWilcox.test(positive ~ test,data=tt.dat,paired=TRUE)


#--------------------------------------
# Spearman's rank correlation coefficient
#--------------------------------------
# CHANGE TO DNA CONCENTRATION
al.corr <- cor.test(x=qdata$CTmean.Al, y=qdata$alepg, 
       method = 'spearman',exact=FALSE)

Na.corr <- cor.test(x=qdata$CTmean.Na, y=qdata$hwkk, 
       method = 'spearman',exact=FALSE)
Ad.corr <- cor.test(x=qdata$CTmean.Ad, y=qdata$hwkk, 
       method = 'spearman',exact=FALSE)
Ac.corr <- cor.test(x=qdata$CTmean.Ac, y=qdata$hwkk, 
       method = 'spearman',exact=FALSE)

tt.corr <- cor.test(x=qdata$CTmean.Tt, y=qdata$ttepg, 
       method = 'spearman',exact=FALSE)

save(al.kk,hw.kk,tt.kk,
     al.q,hw.q,na.q,ac.q,ad.q,tt.q,
     al.kk.gmn,hw.kk.gmn,tt.kk.gmn,
     al.q.gmn,na.q.gmn,ac.q.gmn,ad.q.gmn,tt.q.gmn,
     file="~/Dropbox/WASH-B-STH-Add-on/TFGH/Data/RData/prev_results.RData")
