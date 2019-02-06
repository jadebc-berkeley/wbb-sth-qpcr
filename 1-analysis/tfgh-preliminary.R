#######################################
# WASH Benefits Bangladesh STH KK qPCR validation
# Clean and merge qPCR and KK data
#######################################
rm(list=ls())
library(reshape2)
library(ggplot2)
library(dplyr)
library(tidyr)

#--------------------------------------
# read in qPCR data
#--------------------------------------
data.dir="~/Dropbox/WASH-B-STH-Add-on/TFGH/Data/Final file/"
res.name="FINAL FINAL FINAL Master Results File - Bangladesh - 4-30-18"
iac=read.csv(file=paste0(data.dir,res.name,"-IAC.csv"),stringsAsFactors=FALSE)
na=read.csv(file=paste0(data.dir,res.name,"-NA.csv"),stringsAsFactors=FALSE)
ad=read.csv(file=paste0(data.dir,res.name,"-AD.csv"),stringsAsFactors=FALSE)
ss=read.csv(file=paste0(data.dir,res.name,"-SS.csv"),stringsAsFactors=FALSE)
tt=read.csv(file=paste0(data.dir,res.name,"-TT.csv"),stringsAsFactors=FALSE)
al=read.csv(file=paste0(data.dir,res.name,"-AL.csv"),stringsAsFactors=FALSE)
ac=read.csv(file=paste0(data.dir,res.name,"-AC.csv"),stringsAsFactors=FALSE)

colnames=c("sampleid","sampleno","assay",
  "CTmean","CTSD","assaydate","plateid","comments")

colnames(iac)=colnames
colnames(na)=colnames
colnames(ad)=colnames
colnames(ss)=colnames
colnames(tt)=colnames
colnames(al)=colnames
colnames(ac)=colnames

# drop problematic sample
na=na[na$comments=="",]
ad=ad[ad$comments=="",]
ss=ss[ss$comments=="",]
tt=tt[tt$comments=="",]
ac=ac[ac$comments=="",]

# correcting assay column in Na
na <- na %>% mutate(assay="Na")

qpcr=rbind(iac,na,ad,ss,tt,al,ac)
qpcr=qpcr[!is.na(qpcr$sampleid),]

# check number of results per assay
table(qpcr$assay)

# correcting data entry error
qpcr$sampleid[qpcr$sampleid==":559901ETS"]="559901ETS1"

# take mean of results for each sample
qpcr <- qpcr %>%
  group_by(sampleid,assay) %>%
  summarise(CTmean=mean(CTmean),CTSD=mean(CTSD)) %>%
  mutate(positive=ifelse(CTmean<40,1,0)) %>%
  mutate(positive=ifelse(is.na(CTmean),0,positive))

mean.l <- qpcr %>% 
  group_by(sampleid) %>%
  select(sampleid,assay,CTmean) %>% 
  spread(assay,CTmean) %>%
  rename(CTmean.Ac=Ac, CTmean.Ad=Ad, CTmean.Al=Al,
         CTmean.IAC=IAC, CTmean.Na=Na, 
         CTmean.Ss=Ss,CTmean.Tt=Tt)

sd.l <- qpcr %>% 
  group_by(sampleid) %>%
  select(sampleid,assay,CTSD) %>% 
  spread(assay,CTSD) %>%
  rename(CTSD.Ac=Ac, CTSD.Ad=Ad, CTSD.Al=Al,
         CTSD.IAC=IAC, CTSD.Na=Na, 
         CTSD.Ss=Ss,CTSD.Tt=Tt)

pos.l <- qpcr %>% 
  group_by(sampleid) %>%
  select(sampleid,assay,positive) %>% 
  spread(assay,positive) %>%
  rename(positive.Ac=Ac, positive.Ad=Ad, positive.Al=Al,
         positive.IAC=IAC, positive.Na=Na, 
         positive.Ss=Ss,positive.Tt=Tt) 

qpcr.w <- full_join(mean.l, sd.l, by=c("sampleid"))
qpcr.w <- full_join(qpcr.w, pos.l, by=c("sampleid"))

qpcr.w$dataid=substr(qpcr.w$sampleid,1,5)
qpcr.w$personid=paste(substr(qpcr.w$sampleid,7,7),1,sep="")
qpcr.w$qpcr="Done"

qpcr.w <- qpcr.w %>%
  mutate(positive.Hw=case_when(
    positive.Ac==1 | positive.Na==1 | positive.Ad==1 ~ 1,
    positive.Ac==0 & positive.Na==0 & positive.Ad==0 ~ 0,
    TRUE ~ NA_real_
  )) %>%
  # manual correction of id that doesn't match list sent to Smith
  ungroup() %>%
  mutate(sampleid=ifelse(sampleid=="559901ETS1","59901ETS1",sampleid)) %>%
  mutate(sampleid=ifelse(sampleid=="18705EOS1","18705ECS1",sampleid)) %>%
  mutate(personid=substr(sampleid,7,7),
         dataid=substr(sampleid,1,5)) %>%
  mutate(personid=ifelse(personid=="T","T1",personid)) %>%
  mutate(personid=ifelse(personid=="C","C1",personid)) %>%
  mutate(personid=ifelse(personid=="O","O1",personid)) %>%
  mutate(personid=ifelse(personid=="W","T2",personid)) 


#--------------------------------------
# read in revised Ascaris qPCR data
#--------------------------------------
ascaris_new=read.csv("~/Dropbox/WASH-B-STH-Add-on/TFGH/Data/Revised ascaris assay/KK v New assay comparison_10-17-18.csv",stringsAsFactors=FALSE)
colnames(ascaris_new) = c("sampleid", "alepg", "al", "X", "CTmean.Al2", "CTSD.Al2")
ascaris_new$dataid=substr(ascaris_new$sampleid,1,5)
ascaris_new$personid=paste(substr(ascaris_new$sampleid,7,7),1,sep="")
ascaris_new = ascaris_new %>% select(-c(X, alepg, al))

# indicator for positive Al 
ascaris_new <- ascaris_new %>%
  mutate(CTmean.Al2 = ifelse(CTmean.Al2 == "DNA SAMPLE MISSING", NA,CTmean.Al2)) %>%
  mutate(CTmean.Al2 = as.numeric(CTmean.Al2)) %>%
  mutate(positive.Al2=ifelse(CTmean.Al2<40,1,0)) %>%
  mutate(positive.Al2=ifelse(is.na(CTmean.Al2),0,positive.Al2))

#--------------------------------------
# merge in kk data
#--------------------------------------
kk=read.csv("~/Dropbox/WASHB Parasites/Analysis datasets/Jade/sth.csv")
kk$dataid=as.character(kk$dataid)
kk$personid=as.character(kk$personid)
kk=kk[,c("dataid","personid","block","clusterid","tr",
         "sex","dw","aged","agem","agey",
      "alepg","hwepg","ttepg",
       "logalepg","loghwepg","logttepg",
       "al","tt","hw","sth","alint","ttint","hwint")]

colnames(kk)[which(colnames(kk)=="tt")]="ttkk"
colnames(kk)[which(colnames(kk)=="al")]="alkk"
colnames(kk)[which(colnames(kk)=="hw")]="hwkk"

# merge kk and qPCR data
data=full_join(kk,qpcr.w,by=c("dataid","personid"))

# merge on re-run ascaris qPCR data
data=full_join(data,ascaris_new,by=c("dataid","personid"))

data$qpcr[is.na(data$qpcr)]="Not done"
qdata=data[data$qpcr=="Done",]

# missing characteristics for the 11 kids without kk that got qPCR
missing=qdata[is.na(qdata$alkk),c("dataid","personid")]
missing$dataperson=paste0(missing$dataid,missing$personid)

# clean up missing cluster ids
qdata$clusterid[is.na(qdata$clusterid)]=substr(qdata$dataid[is.na(qdata$clusterid)],1,3)

# -------------------------------------------
# create gold standard by pooling together 
# kk and qPCR results
# -------------------------------------------
qdata = qdata %>% 
  mutate(positive.Hw=ifelse(positive.Na==1 | positive.Ac==1 | positive.Ad==1,1,0)) %>%
  mutate(gold.hwpos=ifelse(hwkk==1 | positive.Hw==1,1,0)) %>%
  mutate(gold.ttpos=ifelse(ttkk==1 | positive.Tt==1,1,0)) %>%
  mutate(gold.alpos=ifelse(alkk==1 | positive.Al==1,1,0)) %>%
  mutate(gold.sthpos=ifelse(sth==1 | positive.Al==1| positive.Hw==1 | positive.Tt==1,1,0))

# -------------------------------------------
# create indicators for moderate/heavy kk infection
# -------------------------------------------
qdata <- qdata %>%
  mutate(almh=ifelse(alepg>=5000,1,0),
         hwmh=ifelse(hwepg>=2000,1,0),
         ttmh=ifelse(ttepg>=1000,1,0),
         almh.f=as.factor(ifelse(almh==1,"Moderate-heavy intensity\ninfection","Low intensity\ninfection")),
         hwmh.f=as.factor(ifelse(hwmh==1,"Moderate-heavy intensity","Low intensity")),
         ttmh.f=as.factor(ifelse(ttmh==1,"Moderate-heavy intensity","Low intensity")))

#--------------------------------------
# save data
#--------------------------------------
save(qdata,file="~/Dropbox/WASH-B-STH-Add-on/TFGH/Data/RData/qdata.RData")
write.csv(qdata,file="~/Dropbox/WASH-B-STH-Add-on/TFGH/Data/RData/qdata.csv",row.names=FALSE)
