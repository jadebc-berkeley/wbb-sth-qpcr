#######################################
# WASH Benefits Bangladesh STH KK qPCR validation
# venn diagram of co-infection
#######################################
rm(list=ls())
library(dplyr)
library(VennDiagram)
library(venneuler)

load("~/Dropbox/WASH-B-STH-Add-on/TFGH/Data/RData/qdata.RData")

# ---------------------------------------
# kato katz plot
# ---------------------------------------
kk=qdata[,c("alkk","hwkk","ttkk")]
kk$coinf=rowSums(kk,na.rm=TRUE)

# intersections
alhw=nrow(kk[kk$alkk==1 & kk$hwkk==1 & kk$ttkk==0,])
hwtt=nrow(kk[kk$alkk==0 & kk$hwkk==1 & kk$ttkk==1,])
altt=nrow(kk[kk$alkk==1 & kk$hwkk==0 & kk$ttkk==1,])
alhwtt=nrow(kk[kk$alkk==1 & kk$hwkk==1 & kk$ttkk==1,])

altot=nrow(kk[kk$alkk==1 & kk$hwkk==0 & kk$ttkk==0,])
hwtot=nrow(kk[kk$alkk==0 & kk$hwkk==1 & kk$ttkk==0,])
tttot=nrow(kk[kk$alkk==0 & kk$hwkk==0 & kk$ttkk==1,])

# circles proportional to size 
v <- venneuler(c(Ascaris=altot, Hookworm=hwtot, Trichuris=tttot,
                 "Ascaris&Hookworm"=alhw,"Hookworm&Trichuris"=hwtt,
                 "Trichuris&Ascaris"=altt,"Ascaris&Hookworm&Trichuris"=alhwtt))
v$labels <- c(
  "","",""
)

cb.lightorange="#E69F00"
cb.blue= "#56B4E9"
cb.green="#009E73"
cb.orange="#D55E00"
cb.pink="#CC79A7"
cb.dblue="#005787"
purple="#9E4AED"

pdf(file="~/Dropbox/WASH-B-STH-Add-on/TFGH/Results/wbb-venn-kk.pdf",
    width=6,height=6)
  plot(v, col=c(cb.lightorange, cb.blue, cb.green))
  # titles
  text(x=v$centers[1,1],y=v$centers[1,2]+0.32, labels=expression(paste(italic("A. lumbricoides"))))
  text(x=v$centers[1,1]+0.3,y=v$centers[1,2]-0.28, labels="Hookworm")
  text(x=v$centers[1,1]+0.3,y=v$centers[1,2]+0.22, labels=expression(paste(italic("T. trichiura"))))
  
  # non intersections
  text(x=v$centers[1,1],y=v$centers[1,2], labels=paste(altot))
  text(x=v$centers[2,1]+0.06,y=v$centers[2,2]-0.05, labels=paste(hwtot))
  text(x=v$centers[3,1]+0.06,y=v$centers[3,2]+0.02, labels=paste(tttot))
  
  # intersections
  text(x=v$centers[1,1]+0.21,y=v$centers[1,2]+0.08, labels=paste(altt))
  text(x=v$centers[1,1]+0.18,y=v$centers[1,2]-0.13, labels=paste(alhw))
  text(x=v$centers[1,1]+0.22,y=v$centers[1,2]-0.02, labels=paste(alhwtt))
  text(x=v$centers[1,1]+0.32,y=v$centers[1,2]-0.04, labels=paste(hwtt))
dev.off()

# ---------------------------------------
# qpcr plot
# ---------------------------------------
# counts
altot=nrow(qdata[qdata$positive.Al2==1 & qdata$positive.Hw==0 & qdata$positive.Tt==0,])
hwtot=nrow(qdata[qdata$positive.Al==0 & qdata$positive.Hw==1 & qdata$positive.Tt==0,])
tttot=nrow(qdata[qdata$positive.Al==0 & qdata$positive.Hw==0 & qdata$positive.Tt==1,])

# intersections
alhw=nrow(qdata[qdata$positive.Al==1 & qdata$positive.Hw==1 & qdata$positive.Tt==0,])
altt=nrow(qdata[qdata$positive.Al==1 & qdata$positive.Hw==0 & qdata$positive.Tt==1,])
hwtt=nrow(qdata[qdata$positive.Al==0 & qdata$positive.Hw==1 & qdata$positive.Tt==1,])

alhwtt=nrow(qdata[qdata$positive.Al2==1 & qdata$positive.Hw==1 & qdata$positive.Tt==1,])

qpcr <- venneuler(c(Ascaris=altot, Hookworm=hwtot, Trichuris=tttot,
                 "Ascaris&Hookworm"=alhw,"Hookworm&Trichuris"=hwtt,
                 "Trichuris&Ascaris"=altt,"Ascaris&Hookworm&Trichuris"=alhwtt))
qpcr$labels <- c("","","")

centers=qpcr$centers

pdf(file="~/Dropbox/WASH-B-STH-Add-on/TFGH/Results/wbb-venn-qpcr.pdf",
    width=6,height=6)
  plot(qpcr, col=c(cb.lightorange, cb.blue,cb.green))
  # titles
  text(x=centers[1,1]+0.15,y=centers[1,2]+0.25, labels=expression(paste(italic("A. lumbricoides"))))
  text(x=centers[2,1]-0.20,y=centers[2,2]+0.25, labels="Hookworm")
  text(x=centers[3,1]-0.19,y=centers[3,2]-0.19, labels=expression(paste(italic("T. trichiura"))))
  
  # non intersections
  text(x=centers[1,1]+0.10,y=centers[1,2]-0.02, labels=paste(altot))
  text(x=centers[2,1]-0.01,y=centers[2,2]+0.05, labels=paste(hwtot))
  text(x=centers[3,1]-0.02,y=centers[3,2]-0.09, labels=paste(tttot))
  
  # intersections
  text(x=centers[1,1]-0.10,y=centers[1,2]-0.12, labels=paste(altt))
  text(x=centers[1,1]-0.16,y=centers[1,2]-0.02, labels=paste(alhwtt))
  text(x=centers[3,1]-0.07,y=centers[3,2]+0.07, labels=paste(hwtt))
  text(x=centers[1,1]-0.1,y=centers[1,2]+0.10, labels=paste(alhw))
dev.off()


# ---------------------------------------
# qpcr hw plot
# ---------------------------------------
# counts
actot=nrow(qdata[qdata$positive.Ac==1 & qdata$positive.Na==0 & qdata$positive.Ad==0,])
natot=nrow(qdata[qdata$positive.Ac==0 & qdata$positive.Na==1 & qdata$positive.Ad==0,])
adtot=nrow(qdata[qdata$positive.Ac==0 & qdata$positive.Na==0 & qdata$positive.Ad==1,])

# intersections
acna=nrow(qdata[qdata$positive.Ac==1 & qdata$positive.Na==1 & qdata$positive.Ad==0,])
acad=nrow(qdata[qdata$positive.Ac==1 & qdata$positive.Na==0 & qdata$positive.Ad==1,])
naad=nrow(qdata[qdata$positive.Ac==0 & qdata$positive.Na==1 & qdata$positive.Ad==1,])
acnaad=nrow(qdata[qdata$positive.Ac==1 & qdata$positive.Na==1 & qdata$positive.Ad==1,])

qpcr.hw <- venneuler(c("A. ceylanicum"=actot, 
                    "N. americanus"=natot, 
                    "A. duodenale"=adtot,
                    "A. ceylanicum&N. americanus"=acna,
                    "A. ceylanicum&A. duodenale"=acad,
                    "N. americanus&A. duodenale"=naad,
                    "A. ceylanicum&N. americanus&A. duodenale"=acnaad))

qpcr.hw$labels <- c("","","")

centers=qpcr.hw$centers

pdf(file="~/Dropbox/WASH-B-STH-Add-on/TFGH/Results/wbb-venn-qpcr-hw.pdf",
    width=6,height=6)
plot(qpcr.hw, col=c(cb.lightorange, cb.blue, purple,cb.green))
# titles
text(x=centers[1,1],y=centers[1,2]+0.15, labels=expression(paste(italic("A. ceylanicum"))))
text(x=centers[2,1],y=centers[2,2]-0.32, labels=expression(paste(italic("N. americanus"))))
text(x=centers[3,1]-0.08,y=centers[3,2]+0.04, labels=expression(paste(italic("A. duodenale"))))
# non intersections
text(x=centers[1,1],y=centers[1,2]+0.03, labels=paste(actot))
text(x=centers[2,1],y=centers[2,2]-0.02, labels=paste(natot))
# intersections
text(x=centers[1,1]+0.01,y=centers[1,2]-0.08, labels=paste(acna))
text(x=centers[1,1]-0.2,y=centers[1,2]-0.2, labels=paste(naad),cex=0.75)
dev.off()
