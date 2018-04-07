#######################################
# WASH Benefits Bangladesh STH KK qPCR validation
# sensitivity and specificity analysis
#######################################
rm(list=ls())
library(dplyr)
library(clusrank)
library(tidyr)
load("~/Dropbox/WASH-B-STH-Add-on/TFGH/Data/RData/qdata.RData")

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

