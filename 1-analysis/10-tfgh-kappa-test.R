#######################################
# WASH Benefits Bangladesh STH KK qPCR validation
# kappa test of agreement between 
# KK and qPCR classification
#######################################
rm(list=ls())
library(irr)

load("~/Dropbox/WASH-B-STH-Add-on/TFGH/Data/RData/qdata.RData")

qdata = qdata %>% filter(!is.na(positive.Al2))

al.kappa = kappa2(qdata[,c("alkk","positive.Al2")], "unweighted")
hw.kappa = kappa2(qdata[,c("hwkk","positive.Hw")], "unweighted")
tt.kappa = kappa2(qdata[,c("ttkk","positive.Tt")], "unweighted")

save(al.kappa,  hw.kappa,  tt.kappa,
        file="~/Dropbox/WASH-B-STH-Add-on/TFGH/Data/RData/kappa_test.RData")
