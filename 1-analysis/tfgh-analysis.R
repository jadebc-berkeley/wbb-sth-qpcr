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
qdata.conc <- qdata.conc %>%
  mutate(copies.Al=ifelse(is.na(copies.Al),0,copies.Al),
         copies.Ac=ifelse(is.na(copies.Ac),0,copies.Ac),
         copies.Na=ifelse(is.na(copies.Na),0,copies.Na),
         copies.Tt=ifelse(is.na(copies.Tt),0,copies.Tt)) 

al.corr <- cor.test(x=qdata.conc$copies.Al, y=qdata.conc$alepg, 
       method = 'spearman',exact=FALSE)
Na.corr <- cor.test(x=qdata.conc$copies.Na, y=qdata.conc$hwepg, 
       method = 'spearman',exact=FALSE)
Ac.corr <- cor.test(x=qdata.conc$copies.Ac, y=qdata.conc$hwepg, 
       method = 'spearman',exact=FALSE)
tt.corr <- cor.test(x=qdata.conc$copies.Tt, y=qdata.conc$ttepg, 
       method = 'spearman',exact=FALSE)

# bootstrap to get p-value for spearman rank correlation
# resample clusters
bs.clus=function(y1, y2,clus,B){
  
  dat=as.data.frame(cbind(y1,y2,clus))
  # dat=dat[dat$clus<=15,]
  # clus=dat$clus
  # y1=dat$y1
  # y2=dat$y2
  
  # observed rho on full data
  corr.obs=cor.test(x=dat$y1, y=dat$y2, 
                method = 'spearman',exact=FALSE)$estimate
  
  # get bootstrap distribution
  bdat=matrix(NA,nrow=B,ncol=1)
  
  blocks=unique(clus)
  
  for(i in 1:B){
    # take a random sample of clusters
    samp.b=sample(blocks,size=length(blocks),replace=TRUE)

    # create data frame with sampled clusters
    bootdat <- NULL
    clus1 <- dat[clus %in% samp.b[1],]
    for(j in 2:length(samp.b)){
      # cc <- dat[clus %in% names(samp.blocks.tab[samp.blocks.tab %in% j]),]
      clusj <- dat[clus %in% samp.b[j],]
      if(j==2){
        bootdat <- rbind(clus1, clusj)
      }else{
        bootdat <- rbind(bootdat, clusj)
      }
    
    # correlation test on bootstrapped data
      corr=cor.test(x=bootdat$y1, y=bootdat$y2, 
                    method = 'spearman',exact=FALSE)
      bdat[i,1] <- corr$estimate
    }
    

  }
  
  # p-value 
  pval = sum(bdat >= corr.obs)/B

  return(list(bsdist=bdat, corr.obs=corr.obs, pval=pval))
}

y1=qdata.conc$alepg[!is.na(qdata.conc$alkk)]
y2=qdata.conc$copies.Al[!is.na(qdata.conc$alkk)]
clus=qdata.conc$block[!is.na(qdata.conc$alkk)]
B=10

set.seed(333)
x=bs.clus(y1=qdata.conc$alepg, y2=qdata.conc$copies.Al, clus=qdata.conc$block,
        B=1000)
x$pval
hist(x$bsdist)

corr.tab=data.frame(org=c("A. lumbricoides","N. americanus","A. ceylanicum",
                          "T. trichiura"),
                    rho=c(al.corr$estimate, Na.corr$estimate, Ac.corr$estimate,
                          tt.corr$estimate),
                    p.val=c(al.corr$p.value, Na.corr$p.value, Ac.corr$p.value,
                            tt.corr$p.value))
corr.tab$rho=sprintf("%0.03f",corr.tab$rho)
corr.tab$p.val=formatC(corr.tab$p.val, format = "e", digits = 2)

write.csv(corr.tab,file="~/Dropbox/WASH-B-STH-Add-on/TFGH/Results/corr.table.csv",row.names=FALSE)

save(al.kk,hw.kk,tt.kk,
     al.q,hw.q,na.q,ac.q,ad.q,tt.q,
     al.kk.gmn,hw.kk.gmn,tt.kk.gmn,
     al.q.gmn,na.q.gmn,ac.q.gmn,ad.q.gmn,tt.q.gmn,
     corr.tab,
     file="~/Dropbox/WASH-B-STH-Add-on/TFGH/Data/RData/prev_results.RData")
