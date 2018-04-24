#######################################
# WASH Benefits Bangladesh STH KK qPCR validation
# venn diagram of co-infection
# try the second package here.. just cant install
# on  my computer

# https://rstudio-pubs-static.s3.amazonaws.com/13301_6641d73cfac741a59c0a851feb99e98b.html

#######################################
rm(list=ls())
library(dplyr)
# library(gdata)
# library(gplots)
library(VennDiagram)
load("~/Dropbox/WASH-B-STH-Add-on/TFGH/Data/RData/qdata.RData")

kk=qdata[,c("alkk","hwkk","ttkk")]
kk$coinf=rowSums(kk,na.rm=TRUE)

# ascaris+hookworm
alhw=nrow(kk[kk$alkk==1 & kk$hwkk==1,])
hwtt=nrow(kk[kk$hwkk==1 & kk$ttkk==1,])
altt=nrow(kk[kk$alkk==1 & kk$ttkk==1,])
alhwtt=nrow(kk[kk$alkk==1 & kk$hwkk==1 & kk$ttkk==1,])

altot=nrow(kk[kk$alkk==1,])
hwtot=nrow(kk[kk$hwkk==1,])
tttot=nrow(kk[kk$ttkk==1,])

draw.triple.venn(area1 = altot, area2 = hwtot, area3 = tttot, 
                 n12 = alhw, n23 = hwtt, n13 = altt, n123 = alhwtt, 
                 category = c("Ascaris", "Hookworm", "Trichuris"), 
                 fill = c("skyblue", "pink1", "mediumorchid"),
                 euler.d=TRUE,scaled=TRUE)


# circles not proportional to size 