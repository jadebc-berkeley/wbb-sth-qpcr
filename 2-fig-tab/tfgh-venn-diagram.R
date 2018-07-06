#######################################
# WASH Benefits Bangladesh STH KK qPCR validation
# venn diagram of co-infection
#######################################
rm(list=ls())
library(dplyr)
library(VennDiagram)
library(venneuler)

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

# circles proportional to size 
v <- venneuler(c(Ascaris=altot, Hookworm=hwtot, Trichuris=tttot,
                 "Ascaris&Hookworm"=alhw,"Hookworm&Trichuris"=hwtt,
                 "Trichuris&Ascaris"=altt,"Ascaris&Hookworm&Trichuris"=alhwtt))
v$labels <- c(
  paste("Ascaris\n",altot),
  paste("Hookworm\n",hwtot),
  paste("Trichuris\n",tttot)
)

pdf(file="~/Dropbox/WASH-B-STH-Add-on/TFGH/Results/wbb-qpcr-kk-venn.pdf",
    width=6,height=6)
  plot(v)
  text(x=v$centers[1,1]+0.2,y=v$centers[1,2]+0.06, labels=paste(altt))
  text(x=v$centers[1,1]+0.15,y=v$centers[1,2]-0.1, labels=paste(alhw))
  text(x=v$centers[1,1]+0.22,y=v$centers[1,2]-0.03, labels=paste(alhwtt))
  text(x=v$centers[1,1]+0.31,y=v$centers[1,2]-0.04, labels=paste(hwtt))
dev.off()



