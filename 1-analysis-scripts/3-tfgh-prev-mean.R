#######################################
# WASH Benefits Bangladesh STH KK qPCR validation
# calculate prevalence and geometric mean
# using each diagnostic

# by Jade Benjamin-Chung
# jadebc@berkeley.edu
#######################################
rm(list=ls())

# configure directories, load libraries and base functions
source(paste0(here::here(), "/0-config.R"))

# load data
load(paste0(data_dir,"qdata.RData"))

#--------------------------------------
# Estimate prevalence and 95%CI for each
# test, robust sandwich SEs
#--------------------------------------
al.kk=washb_mean(qdata$alkk, id=qdata$clusterid, print = TRUE)
hw.kk=washb_mean(qdata$hwkk, id=qdata$clusterid, print = TRUE)
tt.kk=washb_mean(qdata$ttkk, id=qdata$clusterid, print = TRUE)
sth.kk=washb_mean(qdata$sth, id=qdata$clusterid, print = TRUE)

al.q=washb_mean(qdata$positive.Al, id=qdata$clusterid, print = TRUE)
hw.q=washb_mean(qdata$positive.Hw, id=qdata$clusterid, print = TRUE)
na.q=washb_mean(qdata$positive.Na, id=qdata$clusterid, print = TRUE)
ac.q=washb_mean(qdata$positive.Ac, id=qdata$clusterid, print = TRUE)
ad.q=washb_mean(qdata$positive.Ad, id=qdata$clusterid, print = TRUE)
tt.q=washb_mean(qdata$positive.Tt, id=qdata$clusterid, print = TRUE)
ss.q=washb_mean(qdata$positive.Ss, id=qdata$clusterid, print = TRUE)
sth.q=washb_mean(qdata$positive.Sth, id=qdata$clusterid, print = TRUE)

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
prop.table(table(qdata$almh))*100
prop.table(table(qdata$hwmh))*100
prop.table(table(qdata$ttmh))*100

quantile(qdata$CTmean.Al[qdata$almh==1],probs=c(0,0.5,1),na.rm=TRUE)
quantile(qdata$CTmean.Al[qdata$almh==0],probs=c(0,0.5,1),na.rm=TRUE)

#--------------------------------------
# moderate intensity infection
#--------------------------------------
prop.table(table(qdata$al_mod))*100
prop.table(table(qdata$hw_mod))*100
prop.table(table(qdata$tt_mod))*100

quantile(qdata$CTmean.Al[qdata$al_mod==1],probs=c(0,0.5,1),na.rm=TRUE)
quantile(qdata$CTmean.Al[qdata$al_mod==0],probs=c(0,0.5,1),na.rm=TRUE)

#--------------------------------------
# heavy intensity infection
#--------------------------------------
prop.table(table(qdata$al_heavy))*100
prop.table(table(qdata$hw_heavy))*100
prop.table(table(qdata$tt_heavy))*100

quantile(qdata$CTmean.Al[qdata$al_heavy==1],probs=c(0,0.5,1),na.rm=TRUE)
quantile(qdata$CTmean.Al[qdata$al_heavy==0],probs=c(0,0.5,1),na.rm=TRUE)


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

al.kk.gmn.pos=washb_mean(qdata$ln.alepg[qdata$alepg>0], 
                         id=qdata$clusterid[qdata$alepg>0], print = TRUE)
hw.kk.gmn.pos=washb_mean(qdata$ln.hwepg[qdata$hwepg>0], 
                         id=qdata$clusterid[qdata$hwepg>0], print = TRUE)
tt.kk.gmn.pos=washb_mean(qdata$ln.ttepg[qdata$ttepg>0], 
                         id=qdata$clusterid[qdata$ttepg>0], print = TRUE)


al.q.gmn=washb_mean(qdata$ln.CT.al, id=qdata$clusterid, print = TRUE)
na.q.gmn=washb_mean(qdata$ln.CT.na, id=qdata$clusterid, print = TRUE)
ac.q.gmn=washb_mean(qdata$ln.CT.ac, id=qdata$clusterid, print = TRUE)
ad.q.gmn=washb_mean(qdata$ln.CT.ad, id=qdata$clusterid, print = TRUE)
tt.q.gmn=washb_mean(qdata$ln.CT.tt, id=qdata$clusterid, print = TRUE)

#--------------------------------
# exponentiate mean, lb, ub
#--------------------------------
exp_mean_bounds = function(x){
  x = as.data.frame(x)
  x = x %>%
    mutate(Mean = exp(Mean),
           `Lower 95%CI` = exp(`Lower 95%CI`),
           `Upper 95%CI` = exp(`Upper 95%CI`)) %>%
    select(Mean, `Lower 95%CI`, `Upper 95%CI`)
  x =x - 1
  return(x)
}

al.kk.gmn = exp_mean_bounds(al.kk.gmn)
hw.kk.gmn = exp_mean_bounds(hw.kk.gmn)
tt.kk.gmn = exp_mean_bounds(tt.kk.gmn)

al.kk.gmn.pos = exp_mean_bounds(al.kk.gmn.pos)
hw.kk.gmn.pos = exp_mean_bounds(hw.kk.gmn.pos)
tt.kk.gmn.pos = exp_mean_bounds(tt.kk.gmn.pos)

al.q.gmn = exp_mean_bounds(al.q.gmn)
na.q.gmn = exp_mean_bounds(na.q.gmn)
ac.q.gmn = exp_mean_bounds(ac.q.gmn)
tt.q.gmn = exp_mean_bounds(tt.q.gmn)


#--------------------------------------
# Median in positive samples and range
#--------------------------------------
al.med.pos = quantile(qdata$alepg[qdata$alkk==1], probs=c(0,0.5,1))
hw.med.pos = quantile(qdata$hwepg[qdata$hwkk==1], probs=c(0,0.5,1))
tt.med.pos = quantile(qdata$ttepg[qdata$ttkk==1], probs=c(0,0.5,1))

#--------------------------------
# save output
#--------------------------------
save(al.kk,hw.kk,tt.kk,
     al.q,hw.q,na.q,ac.q,ad.q,tt.q,ss.q,
     al.kk.gmn,hw.kk.gmn,tt.kk.gmn,
     al.q.gmn,na.q.gmn,ac.q.gmn,ad.q.gmn,tt.q.gmn,
     al.med.pos, hw.med.pos, tt.med.pos,
     file=paste0(data_dir, "prev_results.RData"))
