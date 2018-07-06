#######################################
# WASH Benefits Bangladesh STH KK qPCR validation
# venn diagram of co-infection
#######################################
rm(list=ls())
library(dplyr)
library(VennDiagram)
library(venneuler)

load("~/Dropbox/WASH-B-STH-Add-on/TFGH/Data/RData/qdata.RData")
load("~/Dropbox/WASH-B-STH-Add-on/TFGH/Data/RData/concentration.RData")

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
  paste("Ascaris\n",altot),
  paste("Hookworm\n",hwtot),
  paste("Trichuris\n",tttot)
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
  text(x=v$centers[1,1]-0.2,y=v$centers[1,2]+0.1, labels=paste(altt))
  text(x=v$centers[1,1]-0.18,y=v$centers[1,2]-0.1, labels=paste(alhw))
  text(x=v$centers[1,1]-0.22,y=v$centers[1,2], labels=paste(alhwtt))
  text(x=v$centers[1,1]-0.34,y=v$centers[1,2]-0.01, labels=paste(hwtt))
dev.off()

# ---------------------------------------
# qpcr plot
# ---------------------------------------
qpcr = qdata.conc %>%
  mutate(positive.Al=ifelse(!is.na(copies.Al),1,0),
         positive.Tt=ifelse(!is.na(copies.Tt),1,0),
         positive.Ss=ifelse(!is.na(copies.Ss),1,0)) %>%
  mutate(positive.Hw=ifelse(positive.Na==1 |positive.Ac==1 | positive.Ad==1,1,0)) %>%
  select(c("positive.Al","positive.Tt","positive.Ss",
           "positive.Hw")) 

# counts
altot=nrow(qpcr[qpcr$positive.Al==1 & qpcr$positive.Hw==0 & qpcr$positive.Tt==0 & qpcr$positive.Ss==0,])
hwtot=nrow(qpcr[qpcr$positive.Al==0 & qpcr$positive.Hw==1 & qpcr$positive.Tt==0 & qpcr$positive.Ss==0,])
tttot=nrow(qpcr[qpcr$positive.Al==0 & qpcr$positive.Hw==0 & qpcr$positive.Tt==1 & qpcr$positive.Ss==0,])
sstot=nrow(qpcr[qpcr$positive.Al==0 & qpcr$positive.Hw==0 & qpcr$positive.Tt==0 & qpcr$positive.Ss==1,])

# intersections
alhw=nrow(qpcr[qpcr$positive.Al==1 & qpcr$positive.Hw==1 & qpcr$positive.Tt==0 & qpcr$positive.Ss==0,])
altt=nrow(qpcr[qpcr$positive.Al==1 & qpcr$positive.Hw==0 & qpcr$positive.Tt==1 & qpcr$positive.Ss==0,])
alss=nrow(qpcr[qpcr$positive.Al==1 & qpcr$positive.Hw==0 & qpcr$positive.Tt==0 & qpcr$positive.Ss==1,])
hwtt=nrow(qpcr[qpcr$positive.Al==0 & qpcr$positive.Hw==1 & qpcr$positive.Tt==1 & qpcr$positive.Ss==0,])
hwss=nrow(qpcr[qpcr$positive.Al==0 & qpcr$positive.Hw==1 & qpcr$positive.Tt==0 & qpcr$positive.Ss==1,])
ttss=nrow(qpcr[qpcr$positive.Al==0 & qpcr$positive.Hw==0 & qpcr$positive.Tt==1 & qpcr$positive.Ss==1,])

alhwtt=nrow(qpcr[qpcr$positive.Al==1 & qpcr$positive.Hw==1 & qpcr$positive.Tt==1 & qpcr$positive.Ss==0,])
alttss=nrow(qpcr[qpcr$positive.Al==1 & qpcr$positive.Hw==0 & qpcr$positive.Tt==1 & qpcr$positive.Ss==1,])
alhwss=nrow(qpcr[qpcr$positive.Al==1 & qpcr$positive.Hw==1 & qpcr$positive.Tt==0 & qpcr$positive.Ss==1,])
hwttss=nrow(qpcr[qpcr$positive.Al==0 & qpcr$positive.Hw==1 & qpcr$positive.Tt==1 & qpcr$positive.Ss==1,])

alhwttss=nrow(qpcr[qpcr$positive.Al==1 & qpcr$positive.Hw==1 & qpcr$positive.Tt==1 & qpcr$positive.Ss==1,])

dna <- venneuler(qpcr)

dna$labels <- c(
  paste("A. lumbricoides\n",altot),
  paste("Hookworm\n",hwtot),
  paste("S. stercoralis\n",sstot),
  paste("T. trichiura\n",tttot)
)

centers=dna$centers

pdf(file="~/Dropbox/WASH-B-STH-Add-on/TFGH/Results/wbb-venn-qpcr.pdf",
    width=6,height=6)
  plot(dna, col=c(cb.lightorange, cb.blue, purple,cb.green))
  text(x=centers[1,1]-0.12,y=centers[1,2]+0.07, labels=paste(altt))
  text(x=centers[4,1]+0.1,y=centers[4,2]+0.01, labels=paste(alhwtt))
  text(x=centers[4,1]-0.02,y=centers[4,2]-0.1, labels=paste(hwtt))
  text(x=centers[2,1]+0.13,y=centers[2,2]+0.09, labels=paste(alhw))
dev.off()


# ---------------------------------------
# qpcr hw plot
# ---------------------------------------
qpcr.hw = qdata.conc %>%
  mutate(positive.Ac=ifelse(!is.na(copies.Ac),1,0),
         positive.Na=ifelse(!is.na(copies.Na),1,0),
         positive.Ad=ifelse(!is.na(copies.Ad),1,0)) %>%
  select(c("positive.Ac","positive.Na","positive.Ad")) 

# counts
actot=nrow(qpcr.hw[qpcr.hw$positive.Ac==1 & qpcr.hw$positive.Na==0 & qpcr.hw$positive.Ad==0,])
natot=nrow(qpcr.hw[qpcr.hw$positive.Ac==0 & qpcr.hw$positive.Na==1 & qpcr.hw$positive.Ad==0,])
adtot=nrow(qpcr.hw[qpcr.hw$positive.Ac==0 & qpcr.hw$positive.Na==0 & qpcr.hw$positive.Ad==1,])

# intersections
acna=nrow(qpcr.hw[qpcr.hw$positive.Ac==1 & qpcr.hw$positive.Na==1 & qpcr.hw$positive.Ad==0,])
acad=nrow(qpcr.hw[qpcr.hw$positive.Ac==1 & qpcr.hw$positive.Na==0 & qpcr.hw$positive.Ad==1,])
naad=nrow(qpcr.hw[qpcr.hw$positive.Ac==0 & qpcr.hw$positive.Na==1 & qpcr.hw$positive.Ad==1,])
acnaad=nrow(qpcr.hw[qpcr.hw$positive.Ac==1 & qpcr.hw$positive.Na==1 & qpcr.hw$positive.Ad==1,])

dna.hw <- venneuler(qpcr.hw)

dna.hw$labels <- c(
  paste("A. ceylanicum\n",actot),
  paste("N. americanus\n",natot),
  paste("A. duodenale\n")
)

centers=dna.hw$centers

pdf(file="~/Dropbox/WASH-B-STH-Add-on/TFGH/Results/wbb-venn-qpcr-hw.pdf",
    width=6,height=6)
plot(dna.hw, col=c(cb.lightorange, cb.blue, purple,cb.green))
text(x=centers[1,1],y=centers[1,2]-0.08, labels=paste(acna))
text(x=centers[1,1]-0.25,y=centers[1,2]-0.16, labels=paste(naad))
dev.off()
