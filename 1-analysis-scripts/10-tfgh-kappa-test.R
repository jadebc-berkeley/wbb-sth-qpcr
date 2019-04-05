#######################################
# WASH Benefits Bangladesh STH KK qPCR validation
# kappa test of agreement between 
# KK and qPCR classification

# by Jade Benjamin-Chung
# jadebc@berkeley.edu
#######################################
rm(list=ls())

# configure directories, load libraries and base functions
source(paste0(here::here(), "/0-config.R"))

# load data
load(paste0(data_dir,"qdata.RData"))

qdata = qdata %>% filter(!is.na(positive.Al2))

al.kappa = kappa2(qdata[,c("alkk","positive.Al2")], "unweighted")
hw.kappa = kappa2(qdata[,c("hwkk","positive.Hw")], "unweighted")
tt.kappa = kappa2(qdata[,c("ttkk","positive.Tt")], "unweighted")

save(al.kappa,  hw.kappa,  tt.kappa,
        file=paste0(data_dir, "kappa_test.RData"))
