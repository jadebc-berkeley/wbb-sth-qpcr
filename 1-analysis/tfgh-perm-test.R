#######################################
# WASH Benefits Bangladesh STH KK qPCR validation
# clustered permutation test
# to assess whether prevalence is different
# between kk and qPCR
#######################################
rm(list=ls())
library(dplyr)
library(clusrank)
library(tidyr)
library(washb)
load("~/Dropbox/WASH-B-STH-Add-on/TFGH/Data/RData/qdata.RData")

al.data <- qdata %>%
  # subset to relevant columns
  select(c(dataid,personid,block,clusterid,positive.Al,alkk)) %>%
  # convert to long format
  gather(test,positive,positive.Al:alkk) %>%
  mutate(test=case_when(
    test=="positive.Al" ~ "qpcr",
    test=="alkk"~ "kk"
  ))

hw.data <- qdata %>%
  # subset to relevant columns
  select(c(dataid,personid,block,clusterid,positive.Hw,hwkk)) %>%
  # convert to long format
  gather(test,positive,positive.Hw:hwkk) %>%
  mutate(test=case_when(
    test=="positive.Hw" ~ "qpcr",
    test=="hwkk"~ "kk"
  ))

tt.data <- qdata %>%
  # subset to relevant columns
  select(c(dataid,personid,block,clusterid,positive.Tt,ttkk)) %>%
  # convert to long format
  gather(test,positive,positive.Tt:ttkk) %>%
  mutate(test=case_when(
    test=="positive.Tt" ~ "qpcr",
    test=="ttkk"~ "kk"
  ))
  
# permutation tests with block as cluster
set.seed(242524)
p.al <- washb_permute(Y=al.data$positive,tr=al.data$test,pair=al.data$block,contrast=c("kk","qpcr"),nreps=100000)

set.seed(242524)
p.hw <- washb_permute(Y=hw.data$positive,tr=hw.data$test,pair=hw.data$block,contrast=c("kk","qpcr"),nreps=100000)

set.seed(242524)
p.tt <- washb_permute(Y=tt.data$positive,tr=tt.data$test,pair=tt.data$block,contrast=c("kk","qpcr"),nreps=100000)


