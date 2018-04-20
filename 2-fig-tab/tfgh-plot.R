#######################################
# WASH Benefits Bangladesh STH KK qPCR validation
# plot qpcr vs. kk data
#######################################
rm(list=ls())
library(dplyr)
library(tidyr)
library(reshape2)
library(ggplot2)
library(gridExtra)
load("~/Dropbox/WASH-B-STH-Add-on/TFGH/Data/RData/qdata.RData")

#--------------------------------------
# bar graph of qPCR CT and KK EPG results
#--------------------------------------
# calculate percent positive
per.pos = qdata %>% 
  # subset columns
  select(alkk,hwkk,ttkk,
         positive.Al,positive.Tt,positive.Hw) %>%
  # calculate percent positive
  summarise_all(funs(mean(.,na.rm=TRUE))) %>%
  gather(lab,per.pos)


# calculate percent positive
n.pos = qdata %>% 
  # subset columns
  select(alkk,hwkk,ttkk,
         positive.Al,positive.Tt,positive.Na,positive.Ad,
         positive.Ac) %>%
  # calculate sum of positives
  summarise_all(funs(sum(.,na.rm=TRUE)))  %>%
  gather(lab,n)
  
bar.data =   
  # merge percent positive
  full_join(n.pos,per.pos,by="lab") %>%
  # add label
  mutate(org=c("Ascaris lumbricoides","Hookworm","Trichuris trichiura",
               "Ascaris lumbricoides","Trichuris trichiura",
               "Necator americanus","Ancylostoma duodenale",
               "Ancylostoma ceylanicum","Hookworm")) %>%
  # sort
  arrange(org) %>%
  # add label
  mutate(orgcat=c(rep("Hookworm",2),rep("Ascaris",2),rep("Hookworm",3),
                  rep("Trichuris",2))) %>%
  mutate(test=c(rep("qPCR",2),rep(c("Kato-Katz","qPCR"),2),
                "qPCR","Kato-Katz","qPCR")) %>%
  mutate(per.f=paste0(sprintf("%0.0f",per.pos*100),"%"))

# manually reassigning hookworm percentage
bar.data$per.f[bar.data$lab=="positive.Na"]="21%"
bar.data$per.f[bar.data$per.f=="NA%"]=NA
bar.data=bar.data[bar.data$lab!="positive.Hw",]
bar.data$org=factor(bar.data$org,levels=c("Ascaris lumbricoides","Hookworm","Ancylostoma ceylanicum",
    "Ancylostoma duodenale","Necator americanus","Trichuris trichiura"))

teal1="#a1dab4"
teal2="#41b6c4"
blue="#2c7fb8"
dblue="#253494"
green="#31a354"
magenta="#c51b8a"

mycol=c(magenta,teal1,teal2,blue,dblue,green)

pdf(file="~/Dropbox/WASH-B-STH-Add-on/TFGH/Results/wbb-qpcr-kk-bargraph.pdf",
    width=10,height=4)
ggplot(bar.data,aes(x=test,y=n,fill=org),col="black")+
  geom_bar(aes(fill=org),stat="identity",colour="black",
           position='stack')+facet_grid(~orgcat)+
  scale_fill_manual("Organism",values=mycol) +
  scale_y_continuous(limits=c(0,1100))+
  geom_text(aes(label=per.f,vjust=c(-0.3,-0.3,-0.3,-0.3,-0.3,-2.3,-0.3,-0.3)))+
  theme_bw()+ylab("Number infected")+xlab("Diagnostic")
dev.off()


#--------------------------------------
# scatter plot of qPCR CT and KK EPG results
#--------------------------------------
# ADD CODE THAT CONVERTS CT VALUE TO DNA CONCENTRATION

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


# al plot
# fix detection limits
qdata$log10alepg=log10(qdata$alepg)
qdata$log10CTmean.Al[is.na(qdata$log10CTmean.Al)]=1

# IS THE Y AXIS THE SAME AS IN EASTON PAPER?
al.plot=ggplot(qdata[qdata$log10CTmean.Al>1 & qdata$log10alepg>1,],
               aes(x=log10alepg,y=log10CTmean.Al))+
  ylab(expression("log"[10]*"DNA concentration (ng/"*mu*"L)"))+
  xlab(expression("log"[10]*italic(" A. lumbricoides")*" mean EPG"))+
  geom_point(alpha=0.5)+
  theme_bw()

# hw plot
qpcr.l=melt(qdata,id.vars=c("dataid","personid"),
            measure.vars=c("CTmean.Na","CTmean.Ad","CTmean.Ac"))
kk.l=melt(qdata,id.vars=c("dataid","personid"),
          measure.vars=c("hwepg"))

qdata.l=merge(qpcr.l,kk.l,by=c("dataid","personid"))
colnames(qdata.l)=c("dataid","personid","qpcr","qpcr.val","kk","kk.val")
qdata.l=qdata.l[qdata.l$kk.val>0,]

hw.plot=ggplot(qdata.l,aes(x=kk.val,y=qpcr.val,group=qpcr))+
  geom_point(aes(color=qpcr),alpha=.5)+
  ylab(expression("log"[10]*"DNA concentration (ng/"*mu*"L)"))+
  xlab(expression("log"[10]*italic(" Hookworm")*" mean EPG"))+
  theme_bw()

# tt plot
# fix detection limits
qdata$log10ttepg=log10(qdata$ttepg)
qdata$log10CTmean.Tt[is.na(qdata$log10CTmean.Tt)]=1

tt.plot=ggplot(qdata[qdata$log10CTmean.Tt>1 & qdata$log10ttepg>1,],
               aes(x=log10ttepg,y=log10CTmean.Tt))+
  geom_point(alpha=0.5)+
  ylab(expression("log"[10]*"DNA concentration (ng/"*mu*"L)"))+
  xlab(expression("log"[10]*italic(" T. trichiura")*" mean EPG"))+
  theme_bw()

ggplot(qdata[qdata$CTmean.Tt>1 & qdata$ttepg>1,],
       aes(x=ttepg,y=CTmean.Tt))+
  geom_point(alpha=0.5)+
  ylab(expression("log"[10]*"DNA concentration (ng/"*mu*"L)"))+
  xlab(expression("log"[10]*italic(" T. trichiura")*" mean EPG"))+
  theme_bw()

# MAKE AXIS LABELS LIKE EASTON'S
cont.plot=grid.arrange(al.plot,hw.plot,tt.plot,nrow=1)
grid.draw(cont.plot)

