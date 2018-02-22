#######################################
# WASH Benefit STH KK qPCR validation
# preliminary analysis
# November 2017
#######################################
rm(list=ls())
library(reshape2)
library(ggplot2)

#--------------------------------------
# read in qPCR data
#--------------------------------------
iac=read.csv(file="~/Dropbox/WASH-B-STH-Add-on/TFGH/Data/File 3/Bangladesh - STH - Master Results File-IAC.csv",stringsAsFactors=FALSE)
na=read.csv(file="~/Dropbox/WASH-B-STH-Add-on/TFGH/Data/File 3/Bangladesh - STH - Master Results File-NA.csv",stringsAsFactors=FALSE)
ad=read.csv(file="~/Dropbox/WASH-B-STH-Add-on/TFGH/Data/File 3/Bangladesh - STH - Master Results File-AD.csv",stringsAsFactors=FALSE)
ss=read.csv(file="~/Dropbox/WASH-B-STH-Add-on/TFGH/Data/File 3/Bangladesh - STH - Master Results File-SS.csv",stringsAsFactors=FALSE)
tt=read.csv(file="~/Dropbox/WASH-B-STH-Add-on/TFGH/Data/File 3/Bangladesh - STH - Master Results File-TT.csv",stringsAsFactors=FALSE)
al=read.csv(file="~/Dropbox/WASH-B-STH-Add-on/TFGH/Data/File 3/Bangladesh - STH - Master Results File-AL.csv",stringsAsFactors=FALSE)
ac=read.csv(file="~/Dropbox/WASH-B-STH-Add-on/TFGH/Data/File 3/Bangladesh - STH - Master Results File-AC.csv",stringsAsFactors=FALSE)

colnames=c("sampleid","sampleno","assay",
  "CTmean","CTSD","assaydate","plateid","retest")

colnames(iac)=colnames
colnames(na)=colnames
colnames(ad)=colnames
colnames(ss)=colnames
colnames(tt)=colnames
colnames(al)=colnames
colnames(ac)=c(colnames,"comment")

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

ac$comment=NULL

qpcr=rbind(iac,na,ad,ss,tt,al,ac)
qpcr=qpcr[!is.na(qpcr$sampleid),]

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



#--------------------------------------
# tabulate data with qPCR and KK results
#--------------------------------------
# qdata$al=ifelse(qdata$al==1,"Yes - qPCR","No - qPCR")
# qdata$alkk=ifelse(qdata$alkk==1,"Yes - KK","No - KK")
# 
# qdata$hw=ifelse(qdata$na==1 | qdata$ad==1 | qdata$ac==1,
#                 "Yes - qPCR","No - qPCR")
# qdata$hwkk=ifelse(qdata$hwkk==1,"Yes - KK","No - KK")
# 
# qdata$tt=ifelse(qdata$tt==1,"Yes - qPCR","No - qPCR")
# # qdata$ttkk=ifelse(qdata$ttkk==1,"Yes - KK","No - KK")
# 
# table(qdata$al,qdata$alkk)
# prop.table(table(qdata$al,qdata$alkk),2)
# 
# # table(qdata$hw,qdata$hwkk)
# # prop.table(table(qdata$tt,qdata$ttkk),2)
# 
# table(qdata$tt,qdata$ttkk)
# prop.table(table(qdata$tt,qdata$ttkk),2)

#--------------------------------------
# bar graph of qPCR CT and KK EPG results
#--------------------------------------
qdata.tot=aggregate(qdata[,c("alkk","hwkk","ttkk",
   "positive.Al","positive.Tt","positive.Na","positive.Ad",
   "positive.Ac")],by=list(qdata$dataid,qdata$personid)
   ,sum)

colnames(qdata.tot)[1:2]=c("dataid","personid")

qdata.l=melt(qdata.tot,id.vars=c("dataid","personid"),
               measure.vars=c("alkk","hwkk","ttkk",
     "positive.Al","positive.Tt","positive.Na","positive.Ad",
                              "positive.Ac"))

n.pos=colSums(qdata[,c("positive.Al","alkk","positive.Na",
      "positive.Ac","positive.Ad","hwkk","positive.Tt","ttkk")],na.rm=TRUE)
N.pos=apply(qdata[,c("positive.Al","alkk","positive.Na",
                     "positive.Ac","positive.Ad","hwkk","positive.Tt","ttkk")],
            2, function(x) length(which(!is.na(x))))
per.pos=n.pos/N.pos

bar.data=data.frame(percent=per.pos,N=N.pos,n=n.pos)
bar.data$org=c("Ascaris lumbricoides","Ascaris lumbricoides",
               "Necator americanus","Ancylostoma ceylanicum",
               "Ancylostoma duodenale","Hookworm",
               "Trichuris trichiura","Trichuris trichiura")
bar.data$orgcat=c("Ascaris","Ascaris",rep("Hookworm",4),
                 rep("Trichuris",2))
bar.data$test=c("qPCR","Kato-Katz","qPCR","qPCR","qPCR",
                "Kato-Katz","qPCR","Kato-Katz")
bar.data$n.f=bar.data$n
bar.data$n.f[bar.data$n==0]=NA

mycol=c("blue","red","blue","red","pink","blue","red")

pdf(file="~/Dropbox/WASH-B-STH-Add-on/Results/wbb-qpcr-prelim.pdf",
    width=10,height=4)
ggplot(bar.data,aes(x=test,y=n,fill=org),col="black")+
  geom_bar(aes(fill=org),stat="identity",
           position='stack')+facet_grid(~orgcat)+
  geom_text(aes(label=n.f,vjust=-0.3))+theme_bw()
dev.off()


#--------------------------------------
# scatter plot of qPCR CT and KK EPG results
#--------------------------------------
qdata$CTmean.Al=as.numeric(as.character(qdata$CTmean.Al))
qdata$CTmean.Ad=as.numeric(as.character(qdata$CTmean.Ad))
qdata$CTmean.Tt=as.numeric(as.character(qdata$CTmean.Tt))
qdata$CTmean.Ac=as.numeric(as.character(qdata$CTmean.Ac))
qdata$CTmean.Na=as.numeric(as.character(qdata$CTmean.Na))

qdata$log10CTmean.Al=log10(qdata$CTmean.Al)
qdata$log10CTmean.Ad=log10(qdata$CTmean.Ad)
qdata$log10CTmean.Tt=log10(qdata$CTmean.Tt)
qdata$log10CTmean.Ac=log10(qdata$CTmean.Ac)
qdata$log10CTmean.Na=log10(qdata$CTmean.Na)

# fix detection limits
qdata$log10alepg=log10(qdata$alepg)

qdata$log10CTmean.Al[is.na(qdata$log10CTmean.Al)]=1

# al plot
al.plot=ggplot(qdata,aes(x=log10alepg,y=log10CTmean.Al))+
  geom_point(aes(alpha=0.5))+

# hw plot
qpcr.l=melt(qdata,id.vars=c("dataid","personid"),
        measure.vars=c("CTmean.Na","CTmean.Ad","CTmean.Ac"))
kk.l=melt(qdata,id.vars=c("dataid","personid"),
            measure.vars=c("hwepg"))

qdata.l=merge(qpcr.l,kk.l,by=c("dataid","personid"))
colnames(qdata.l)=c("dataid","personid","qpcr","qpcr.val","kk","kk.val")
qdata.l=qdata.l[qdata.l$kk.val>0,]

hw.plot=ggplot(qdata.l,aes(x=kk.val,y=qpcr.val,group=qpcr))+
  geom_point(aes(color=qpcr,alpha=.5))

# tt plot
tt.plot=ggplot(qdata[qdata$ttepg>0,],aes(x=ttepg,y=CTmean.Tt))+
  geom_point(aes(alpha=0.5))+
  scale_x_continuous(limits=c(0,200))

# select subsample for DNA barcoding / sequencing
# no hookworm or tt
# many ascaris eggs
qdata[(qdata$alint=="High intensity" | 
             qdata$alint=="Moderate intensity") &
            !is.na(qdata$alint),c("dataid","personid","alepg")]
