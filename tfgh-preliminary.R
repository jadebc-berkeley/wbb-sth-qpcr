#######################################
# WASH Benefits Bangladesh STH KK qPCR validation
# Clean data
#######################################
rm(list=ls())
library(reshape2)
library(ggplot2)

#--------------------------------------
# read in qPCR data
#--------------------------------------
data.dir="~/Dropbox/WASH-B-STH-Add-on/TFGH/Data/Full File/"
res.name="Bangladesh - STH - Master Results File - with Retests"
iac=read.csv(file=paste0(data.dir,res.name,"-IAC.csv"),stringsAsFactors=FALSE)
na=read.csv(file=paste0(data.dir,res.name,"-NA.csv"),stringsAsFactors=FALSE)
ad=read.csv(file=paste0(data.dir,res.name,"-AD.csv"),stringsAsFactors=FALSE)
ss=read.csv(file=paste0(data.dir,res.name,"-SS.csv"),stringsAsFactors=FALSE)
tt=read.csv(file=paste0(data.dir,res.name,"-TT.csv"),stringsAsFactors=FALSE)
al=read.csv(file=paste0(data.dir,res.name,"-AL.csv"),stringsAsFactors=FALSE)
ac=read.csv(file=paste0(data.dir,res.name,"-AC.csv"),stringsAsFactors=FALSE)

colnames=c("sampleid","sampleno","assay",
  "CTmean","CTSD","assaydate","plateid","retest")

colnames(iac)=colnames
colnames(na)=colnames
colnames(ad)=colnames
colnames(ss)=colnames
colnames(tt)=colnames
colnames(al)=c(colnames,"comment")
colnames(ac)=c(colnames,"comment")


# WHAT DOES IT MEAN THAT THERE ARE IAC RESULTS IN THE NA
# SPREADSHEET? FOR now assuming they should be NA
na$assay="Na"

# Nils: First off, we generally call any samples that give a positive result with Ct values <40 
# positive if they are positive in both replicate reactions.  If only one replicate is positive, 
# then the sample is retested by another two replicates.  If at least one of these is positive 
# (meaning at least 2 out of 4 overall) the sample is considered to be positive.  
iac$positive=ifelse(iac$CTmean<40 & is.na(iac$retest),1,0)
na$positive=ifelse(na$CTmean<40 & is.na(na$retest),1,0)
ad$positive=ifelse(ad$CTmean<40 & is.na(ad$retest),1,0)
ss$positive=ifelse(ss$CTmean<40 & is.na(ss$retest),1,0)
tt$positive=ifelse(tt$CTmean<40 & is.na(tt$retest),1,0)
al$positive=ifelse(al$CTmean<40 & is.na(al$retest),1,0)
ac$positive=ifelse(ac$CTmean<40 & is.na(ac$retest),1,0)

iac$positive[is.na(iac$positive)]=0 #retest 2 less than theirs
na$positive[is.na(na$positive)]=0
ad$positive[is.na(ad$positive)]=0
ss$positive[is.na(ss$positive)]=0
tt$positive[is.na(tt$positive)]=0
al$positive[is.na(al$positive)]=0 # positive count under theirs by 1
ac$positive[is.na(ac$positive)]=0

al$comment=NULL
ac$comment=NULL

qpcr=rbind(iac,na,ad,ss,tt,al,ac)
qpcr=qpcr[!is.na(qpcr$sampleid),]
qpcr=qpcr[,c("sampleid","assay",
        "CTmean","CTSD","positive")]

qpcr.w=reshape(qpcr[,c("sampleid","assay",
  "CTmean","CTSD","positive")],idvar="sampleid",
  direction="wide",timevar="assay")

qpcr.w=qpcr.w[!is.na(qpcr.w$sampleid),]
qpcr.w$dataid=substr(qpcr.w$sampleid,1,5)
qpcr.w$personid=paste(substr(qpcr.w$sampleid,7,7),1,sep="")
qpcr.w$qpcr="Done"


#--------------------------------------
# merge in kk data
#--------------------------------------
kk=read.csv("~/Dropbox/WASHB Parasites/Analysis datasets/Jade/sth.csv")
kk=kk[,c("dataid","personid","block","clusterid","tr",
      "alepg","hwepg","ttepg",
       "logalepg","loghwepg","logttepg",
       "al","tt","hw","sth","alint","ttint","hwint")]

colnames(kk)[which(colnames(kk)=="tt")]="ttkk"
colnames(kk)[which(colnames(kk)=="al")]="alkk"
colnames(kk)[which(colnames(kk)=="hw")]="hwkk"

data=merge(kk,qpcr.w,by=c("dataid","personid"),all.x=TRUE,all.y=TRUE)
data$qpcr[is.na(data$qpcr)]="Not done"
qdata=data[data$qpcr=="Done",]

#--------------------------------------
# save data
#--------------------------------------
save(qdata,file="~/Dropbox/WASH-B-STH-Add-on/TFGH/Data/RData/qdata.RData")

#--------------------------------------
# ids with positive ascaris kk and 
# all negative for qPCR
#--------------------------------------
qdata$keep=ifelse(qdata$alepg>0 & qdata$positive.Ad==0 &
  qdata$positive.Na==0 & qdata$positive.Tt==0 & 
  qdata$positive.Al==0 & qdata$positive.Ac==0 &
  qdata$hwepg==0 & qdata$ttepg==0,1,0)

keep=qdata[qdata$keep==1,]
keep=keep[!is.na(keep$alint),]
keep=keep[rev(order(keep$alepg)),]
keep[1:6,c("dataid","personid")]
keep[7:16,c("dataid","personid")]
keep[17:29,c("dataid","personid")]

# requested by Nils on 2/22/18
# kk and qpcr positive for tt and nothing else
qdata$keep2=ifelse(qdata$ttkk==1 & qdata$hwkk==0 & qdata$alkk==0 & 
      qdata$positive.Tt==1 & qdata$positive.Na!=1 & qdata$positive.Ad!=1 & 
      qdata$positive.Ss!=1 & qdata$positive.Al!=1 & qdata$positive.Ac!=1,1,0)
as.matrix(paste0(qdata$dataid[qdata$keep2==1],"E",qdata$personid[qdata$keep2==1],"S1"))

# kk and qpcr positive for Na and nothing else
qdata$keep3=ifelse(qdata$ttkk==0 & qdata$hwkk==1 & qdata$alkk==0 & 
       qdata$positive.Tt!=1 & qdata$positive.Na==1 & qdata$positive.Ad!=1 & 
       qdata$positive.Ss!=1 & qdata$positive.Al!=1 & qdata$positive.Ac!=1,1,0)
as.matrix(paste0(qdata$dataid[qdata$keep3==1],"E",qdata$personid[qdata$keep3==1],"S1"))

# kk and qpcr positive for Al and nothing else
qdata$keep4=ifelse(qdata$ttkk==0 & qdata$hwkk==0 & qdata$alkk==1 & 
       qdata$positive.Tt!=1 & qdata$positive.Na!=1 & qdata$positive.Ad!=1 & 
       qdata$positive.Ss!=1 & qdata$positive.Al==1 & qdata$positive.Ac!=1,1,0)
as.matrix(paste0(qdata$dataid[qdata$keep4==1],"E",qdata$personid[qdata$keep4==1],"S1"))

# kk and qpcr negative for all
qdata$keep5=ifelse(qdata$ttkk==0 & qdata$hwkk==0 & qdata$alkk==0 & 
      qdata$positive.Tt!=1 & qdata$positive.Na!=1 & qdata$positive.Ad!=1 & 
      qdata$positive.Ss!=1 & qdata$positive.Al!=1 & qdata$positive.Ac!=1,1,0)
as.matrix(paste0(qdata$dataid[qdata$keep5==1],"E",qdata$personid[qdata$keep5==1],"S1")[1:15])

