#######################################
# WASH Benefits Bangladesh STH KK qPCR validation
# plot qpcr vs. kk data
#######################################
rm(list=ls())
library(reshape2)
library(ggplot2)
load("~/Dropbox/WASH-B-STH-Add-on/TFGH/Data/RData/qdata.RData")

#--------------------------------------
# bar graph of qPCR CT and KK EPG results
#--------------------------------------
qdata.tot=aggregate(qdata[,c("alkk","hwkk","ttkk",
      "positive.Al","positive.Tt","positive.Na","positive.Ad",
      "positive.Ac")],by=list(qdata$dataid,qdata$personid),sum)

qdata.tot = qdata %>% 
  group_by(dataid,personid) %>%
  select(dataid,personid,alkk,hwkk,ttkk,
         positive.Al,positive.Tt,positive.Na,positive.Ad,
         positive.Ac) %>%
  summarise(alkk=sum(alkk),hwkk=sum(hwkk),)

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
# IS THE Y AXIS THE SAME AS IN EASTON PAPER?
al.plot=ggplot(qdata[qdata$log10CTmean.Al>1 & qdata$log10alepg>1,],
               aes(x=log10alepg,y=log10CTmean.Al))+
  geom_point(aes(alpha=0.5))


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
tt.plot=ggplot(qdata,aes(x=ttepg,y=CTmean.Tt))+
  geom_point(aes(alpha=0.5))

# select subsample for DNA barcoding / sequencing
# no hookworm or tt
# many ascaris eggs
qdata[(qdata$alint=="High intensity" | 
         qdata$alint=="Moderate intensity") &
        !is.na(qdata$alint),c("dataid","personid","alepg")]
