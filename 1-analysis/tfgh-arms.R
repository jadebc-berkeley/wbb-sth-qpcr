#######################################
# WASH Benefits Bangladesh STH KK qPCR validation
# between-arm comparisons
#######################################

# USE TMLE INSTEAD?

rm(list=ls())
library(dplyr)
library(washb)

load("~/Dropbox/WASH-B-STH-Add-on/TFGH/Data/RData/qdata.RData")
source("~/documents/crg/wash-benefits/bangladesh/src/sth/analysis/0-base-programs.R")

# pre-process data
d=qdata
d$tr=droplevels(d$tr)

# drop children missing KK for all sth outcomes
d$allsth=ifelse(!is.na(d$positive.Al) & !is.na(d$positive.Na) 
    & !is.na(d$positive.Ac)  & !is.na(d$positive.Ad)  & 
      !is.na(d$positive.Tt),1,0)
d=d[d$allsth==1,]
d$allsth=NULL

# subset to columns needed for unadjusted PR
df = d %>% 
  select(block,clusterid,tr,positive.Al,positive.Al2,positive.Na,positive.Ac,
         positive.Ad,positive.Hw,positive.Tt) %>%
  rename(al=positive.Al,al2=positive.Al2,
         na=positive.Na,ac=positive.Ac,
         ad=positive.Ad,hw=positive.Hw,tt=positive.Tt) %>%
  mutate(block=as.factor(block),
         sth=ifelse(al==1|na==1|ac==1|ad==1|hw==1|tt==1,1,0)) 


#----------------------------------------------
# n and N for prevalence by arm
#----------------------------------------------
N.al=table(d$tr[!is.na(df$al)])
N.hw=table(d$tr[!is.na(df$hw)])
N.na=table(d$tr[!is.na(df$na)])
N.ac=table(d$tr[!is.na(df$ac)])
N.tt=table(d$tr[!is.na(df$tt)])
N.sth=table(d$tr[!is.na(df$sth)])

n.al=table(df$al,df$tr)[2,]
n.hw=table(df$hw,df$tr)[2,]
n.na=table(df$na,df$tr)[2,]
n.ac=table(df$ac,df$tr)[2,]
n.tt=table(df$tt,df$tr)[2,]
n.sth=table(df$sth,df$tr)[2,]

psth_n_prev_j=data.frame(cbind(n.al,N.al,n.hw,N.hw,n.na,N.na,
     n.ac,N.ac,n.tt,N.tt,n.sth,N.sth))

#----------------------------------------------
# prevalence by arm and time point
#----------------------------------------------
al_prev=t(sapply(levels(df$tr), function(x) washb_mean(
  Y=df$al[df$tr==x],id=df$clusterid[df$tr==x],print=FALSE)))
al_prev2=t(sapply(levels(df$tr), function(x) washb_mean(
  Y=df$al2[df$tr==x],id=df$clusterid[df$tr==x],print=FALSE)))
hw_prev=t(sapply(levels(df$tr), function(x) washb_mean(
  Y=df$hw[df$tr==x],id=df$clusterid[df$tr==x],print=FALSE)))
na_prev=t(sapply(levels(df$tr), function(x) washb_mean(
  Y=df$na[df$tr==x],id=df$clusterid[df$tr==x],print=FALSE)))
ac_prev=t(sapply(levels(df$tr), function(x) washb_mean(
  Y=df$ac[df$tr==x],id=df$clusterid[df$tr==x],print=FALSE)))
tt_prev=t(sapply(levels(df$tr), function(x) washb_mean(
  Y=df$tt[df$tr==x],id=df$clusterid[df$tr==x],print=FALSE)))
sth_prev=t(sapply(levels(df$tr), function(x) washb_mean(
  Y=df$sth[df$tr==x],id=df$clusterid[df$tr==x],print=FALSE)))

colnames(al_prev)=c("N","Prev","SD","Robust SE","lb","ub")
colnames(hw_prev)=c("N","Prev","SD","Robust SE","lb","ub")
colnames(na_prev)=c("N","Prev","SD","Robust SE","lb","ub")
colnames(ac_prev)=c("N","Prev","SD","Robust SE","lb","ub")
colnames(tt_prev)=c("N","Prev","SD","Robust SE","lb","ub")
colnames(sth_prev)=c("N","Prev","SD","Robust SE","lb","ub")

psth_prev_j=cbind(al_prev[,c("Prev","lb","ub")],
                  hw_prev[,c("Prev","lb","ub")],
                  na_prev[,c("Prev","lb","ub")],
                  ac_prev[,c("Prev","lb","ub")],
                  tt_prev[,c("Prev","lb","ub")],
                  sth_prev[,c("Prev","lb","ub")])

colnames(psth_prev_j)=c("al-prev","al-lb","al-ub",
                        "hw-prev","hw-lb","hw-ub",
                        "na-prev","na-lb","na-ub",
                        "ac-prev","ac-lb","ac-ub",
                        "tt-prev","tt-lb","tt-ub",
                        "sth-prev","sth-lb","sth-ub")
#----------------------------------------------
# H1: Unadjusted prevalence ratios; each arm vs. 
# control. PR, CI, P-value
#----------------------------------------------
trlist=c("Water","Sanitation", "WSH")

# Poisson regression for RRs
al_rr_h1_unadj_j=t(apply(matrix(trlist), 1,function(x) washb_mh(Y=df$al,tr=df$tr,strat=df$block,
     contrast=c("Control",x),measure="RR")))

al2_rr_h1_unadj_j=t(apply(matrix(trlist), 1,function(x) washb_mh(Y=df$al2,tr=df$tr,strat=df$block,
      contrast=c("Control",x),measure="RR")))

hw_rr_h1_unadj_j=t(apply(matrix(trlist),1 ,function(x) washb_mh(Y=df$hw,tr=df$tr,strat=df$block,
     contrast=c("Control",x),measure="RR")))

na_rr_h1_unadj_j=t(apply(matrix(trlist),1 ,function(x) washb_mh(Y=df$na,tr=df$tr,strat=df$block,
     contrast=c("Control",x),measure="RR")))

ac_rr_h1_unadj_j=t(apply(matrix(trlist),1 ,function(x) washb_mh(Y=df$ac,tr=df$tr,strat=df$block,
     contrast=c("Control",x),measure="RR")))

tt_rr_h1_unadj_j=t(apply(matrix(trlist),1 ,function(x) washb_mh(Y=df$tt,tr=df$tr,strat=df$block,
     contrast=c("Control",x),measure="RR")))

sth_rr_h1_unadj_j=t(apply(matrix(trlist),1 ,function(x) washb_mh(Y=df$sth,tr=df$tr,strat=df$block,
     contrast=c("Control",x),measure="RR")))

# Linear regression for RDs
al_rd_h1_unadj_j=t(apply(matrix(trlist),1 ,function(x) washb_mh(Y=df$al,tr=df$tr,strat=df$block,
     contrast=c("Control",x),measure="RD")))

hw_rd_h1_unadj_j=t(apply(matrix(trlist),1 ,function(x) washb_mh(Y=df$hw,tr=df$tr,strat=df$block,
     contrast=c("Control",x),measure="RD")))

na_rd_h1_unadj_j=t(apply(matrix(trlist),1 ,function(x) washb_mh(Y=df$na,tr=df$tr,strat=df$block,
     contrast=c("Control",x),measure="RD")))

ac_rd_h1_unadj_j=t(apply(matrix(trlist),1 ,function(x) washb_mh(Y=df$ac,tr=df$tr,strat=df$block,
     contrast=c("Control",x),measure="RD")))

tt_rd_h1_unadj_j=t(apply(matrix(trlist),1 ,function(x) washb_mh(Y=df$tt,tr=df$tr,strat=df$block,
     contrast=c("Control",x),measure="RD")))

sth_rd_h1_unadj_j=t(apply(matrix(trlist),1 ,function(x) washb_mh(Y=df$sth,tr=df$tr,strat=df$block,
     contrast=c("Control",x),measure="RD")))

rownames(al_rr_h1_unadj_j)=c("Water vs C","Sanitation vs C", "WSH vs C")
rownames(hw_rr_h1_unadj_j)=c("Water vs C","Sanitation vs C", "WSH vs C")
rownames(na_rr_h1_unadj_j)=c("Water vs C","Sanitation vs C", "WSH vs C")
rownames(ac_rr_h1_unadj_j)=c("Water vs C","Sanitation vs C", "WSH vs C")
rownames(tt_rr_h1_unadj_j)=c("Water vs C","Sanitation vs C", "WSH vs C")
rownames(sth_rr_h1_unadj_j)=c("Water vs C","Sanitation vs C", "WSH vs C")

rownames(al_rd_h1_unadj_j)=c("Water vs C","Sanitation vs C", "WSH vs C")
rownames(hw_rd_h1_unadj_j)=c("Water vs C","Sanitation vs C", "WSH vs C")
rownames(na_rd_h1_unadj_j)=c("Water vs C","Sanitation vs C", "WSH vs C")
rownames(ac_rd_h1_unadj_j)=c("Water vs C","Sanitation vs C", "WSH vs C")
rownames(tt_rd_h1_unadj_j)=c("Water vs C","Sanitation vs C", "WSH vs C")
rownames(sth_rd_h1_unadj_j)=c("Water vs C","Sanitation vs C", "WSH vs C")


#----------------------------------------------
# H2: Unadjusted prevalence ratios; combined WSH vs. 
# single arms.  PR, CI, P-value
#----------------------------------------------

trlist=c("Water","Sanitation")

# Poisson regression for RRs
al_rr_h2_unadj_j=t(apply(matrix(trlist),1 ,function(x) washb_mh(Y=df$al,tr=df$tr,strat=df$block,
     contrast=c(x,"WSH"),measure="RR")))

hw_rr_h2_unadj_j=t(apply(matrix(trlist),1 ,function(x) washb_mh(Y=df$hw,tr=df$tr,strat=df$block,
     contrast=c(x,"WSH"),measure="RR")))

na_rr_h2_unadj_j=t(apply(matrix(trlist),1 ,function(x) washb_mh(Y=df$na,tr=df$tr,strat=df$block,
     contrast=c(x,"WSH"),measure="RR")))

ac_rr_h2_unadj_j=t(apply(matrix(trlist),1 ,function(x) washb_mh(Y=df$ac,tr=df$tr,strat=df$block,
     contrast=c(x,"WSH"),measure="RR")))

tt_rr_h2_unadj_j=t(apply(matrix(trlist),1 ,function(x) washb_mh(Y=df$tt,tr=df$tr,strat=df$block,
     contrast=c(x,"WSH"),measure="RR")))

sth_rr_h2_unadj_j=t(apply(matrix(trlist),1 ,function(x) washb_mh(Y=df$sth,tr=df$tr,strat=df$block,
     contrast=c(x,"WSH"),measure="RR")))

# Linear regression for RDs
al_rd_h2_unadj_j=t(apply(matrix(trlist),1 ,function(x) washb_mh(Y=df$al,tr=df$tr,strat=df$block,
     contrast=c(x,"WSH"),measure="RD")))

hw_rd_h2_unadj_j=t(apply(matrix(trlist),1 ,function(x) washb_mh(Y=df$hw,tr=df$tr,strat=df$block,
     contrast=c(x,"WSH"),measure="RD")))

na_rd_h2_unadj_j=t(apply(matrix(trlist),1 ,function(x) washb_mh(Y=df$na,tr=df$tr,strat=df$block,
     contrast=c(x,"WSH"),measure="RD")))

ac_rd_h2_unadj_j=t(apply(matrix(trlist),1 ,function(x) washb_mh(Y=df$ac,tr=df$tr,strat=df$block,
     contrast=c(x,"WSH"),measure="RD")))

tt_rd_h2_unadj_j=t(apply(matrix(trlist),1 ,function(x) washb_mh(Y=df$tt,tr=df$tr,strat=df$block,
     contrast=c(x,"WSH"),measure="RD")))

sth_rd_h2_unadj_j=t(apply(matrix(trlist),1 ,function(x) washb_mh(Y=df$sth,tr=df$tr,strat=df$block,
     contrast=c(x,"WSH"),measure="RD")))

rownames(al_rr_h2_unadj_j)=c("WSH vs Water","WSH vs Sanitation")
rownames(hw_rr_h2_unadj_j)=c("WSH vs Water","WSH vs Sanitation")
rownames(na_rr_h2_unadj_j)=c("WSH vs Water","WSH vs Sanitation")
rownames(ac_rr_h2_unadj_j)=c("WSH vs Water","WSH vs Sanitation")
rownames(tt_rr_h2_unadj_j)=c("WSH vs Water","WSH vs Sanitation")
rownames(sth_rr_h2_unadj_j)=c("WSH vs Water","WSH vs Sanitation")

rownames(al_rd_h2_unadj_j)=c("WSH vs Water","WSH vs Sanitation")
rownames(hw_rd_h2_unadj_j)=c("WSH vs Water","WSH vs Sanitation")
rownames(na_rr_h2_unadj_j)=c("WSH vs Water","WSH vs Sanitation")
rownames(ac_rr_h2_unadj_j)=c("WSH vs Water","WSH vs Sanitation")
rownames(tt_rd_h2_unadj_j)=c("WSH vs Water","WSH vs Sanitation")
rownames(sth_rd_h2_unadj_j)=c("WSH vs Water","WSH vs Sanitation")


#----------------------------------------------
# save objects
#----------------------------------------------

save(psth_n_prev_j,psth_prev_j,
     al_prev,hw_prev,na_prev,ac_prev,tt_prev,sth_prev,
     
     al_rr_h1_unadj_j,hw_rr_h1_unadj_j,tt_rr_h1_unadj_j,sth_rr_h1_unadj_j,
	   na_rr_h1_unadj_j,ac_rr_h1_unadj_j,
     al_rd_h1_unadj_j,hw_rd_h1_unadj_j,tt_rd_h1_unadj_j,sth_rd_h1_unadj_j,
     na_rd_h1_unadj_j,ac_rd_h1_unadj_j,
     
     al_rr_h2_unadj_j,hw_rr_h2_unadj_j,tt_rr_h2_unadj_j,sth_rr_h2_unadj_j,
     na_rr_h2_unadj_j,ac_rr_h2_unadj_j,
     al_rd_h2_unadj_j,hw_rd_h2_unadj_j,tt_rd_h2_unadj_j,sth_rd_h2_unadj_j,
     na_rd_h2_unadj_j,ac_rd_h2_unadj_j,
     
     file="~/Dropbox/WASH-B-STH-Add-on/TFGH/Data/RData/tfgh_arms.RData")
