#######################################
# WASH Benefits Bangladesh STH KK qPCR validation
# plot qpcr vs. kk data
#######################################
rm(list=ls())
library(dplyr)
library(tidyr)
library(reshape2)
library(ggplot2)
library(grid)
library(gridExtra)
load("~/Dropbox/WASH-B-STH-Add-on/TFGH/Data/RData/qdata.RData")
load("~/Dropbox/WASH-B-STH-Add-on/TFGH/Data/RData/concentration.RData")

#--------------------------------------
# bar graph of qPCR CT and KK EPG results
#--------------------------------------
# calculate percent positive
per.pos = qdata %>% 
  # subset columns
  select(alkk,hwkk,ttkk,
         positive.Al,positive.Tt,positive.Hw,
         positive.Ac,positive.Na,positive.Ad) %>%
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
  mutate(per.pos=per.pos*100,
         per.f=paste0(sprintf("%0.1f",per.pos),"%"))

# manually reassigning hookworm percentage
bar.data$per.f[bar.data$lab=="positive.Na"]=
  paste0(sprintf("%0.1f",mean(qdata$positive.Hw*100,na.rm=TRUE)),"%")
bar.data=bar.data[bar.data$lab!="positive.Hw",]
bar.data$per.f[bar.data$org=="Ancylostoma duodenale"]=""
bar.data$per.f[bar.data$org=="Ancylostoma ceylanicum"]=""
bar.data$org=factor(bar.data$org,levels=c("Ascaris lumbricoides","Hookworm","Ancylostoma ceylanicum",
    "Ancylostoma duodenale","Necator americanus","Trichuris trichiura"))

cb.lightorange="#E69F00"
cb.blue= "#56B4E9"
cb.green="#009E73"
cb.orange="#D55E00"
cb.pink="#CC79A7"
cb.dblue="#005787"
cb.lblue="#A0D9FA"
teal2="#41b6c4"
purple="#B677E6"
gray="#919191"

mycol=c(cb.lightorange,purple,cb.blue,cb.pink,cb.dblue,cb.green)

pdf(file="~/Dropbox/WASH-B-STH-Add-on/TFGH/Results/wbb-qpcr-kk-bargraph.pdf",
    width=10,height=4)
ggplot(bar.data,aes(x=test,y=per.pos,fill=org),col="black")+
  geom_bar(aes(fill=org),stat="identity",colour="black",
           position='stack')+facet_grid(~orgcat)+
  scale_fill_manual("Organism",values=mycol) +
  scale_y_continuous(limits=c(0,40))+
  geom_text(aes(label=per.f,vjust=c(-0.3,-0.3,-0.3,-0.3,-0.3,-2.9,-0.3,-0.3)))+
  theme_bw()+ylab("Prevalence")+xlab("Diagnostic method")
dev.off()


#--------------------------------------
# scatter plot of qPCR CT and KK EPG results
#--------------------------------------
cb.lightorange="#E69F00"
cb.blue= "#56B4E9"
cb.green="#009E73"
cb.orange="#D55E00"
cb.pink="#CC79A7"
cb.dblue="#005787"

yseq=seq(10,55,5)
xseq=c(1, 10, 100, 1000, 10000, 100000)

# al plot
al <- qdata.conc %>% filter(positive.Al==1 | alkk==1) %>%
  # impute 1 for negative values of epg
  mutate(alepg=ifelse(alepg==0, 1, alepg))

al.plot=ggplot(al, aes(x=alepg, y=CTmean.Al))+
  geom_point(alpha=0.65,col=cb.lightorange)+
  geom_smooth(se=FALSE,col="black",data=al[al$alepg>1,], method='loess')+
  scale_y_log10(labels=yseq, breaks=yseq, limits=c(10^1.1, 10^1.7)) +
  scale_x_log10(labels=xseq, breaks=xseq, limits=c(1, 10^5))+
  xlab(expression(paste("log"[10], italic(" A. lumbricoides"), " EPG")))+
  ylab(expression(paste("log"[10], italic(" A. lumbricoides")," Ct value")))+
  theme_bw()+ggtitle(expression(paste(italic("A. lumbricoides"))))+
  theme(plot.title = element_text(hjust = 0.5))

# hw plot
hw <- qdata.conc %>% filter(positive.Hw==1 | hwkk==1) %>%
  select(c(CTmean.Na, CTmean.Ac, CTmean.Ad, hwepg)) %>%
  gather(hw.species,CT,CTmean.Na:CTmean.Ad) %>%
  mutate(Species=case_when(
    hw.species=="CTmean.Ad" ~ "Ancylostoma duodenale",
    hw.species=="CTmean.Ac" ~ "Ancylostoma ceylanicum",
    hw.species=="CTmean.Na" ~ "Necator americanus"
  ))%>%
  # impute 1 for negative values of epg
  mutate(hwepg=ifelse(hwepg==0, 1, hwepg)) %>%
  mutate(Species=factor(Species, levels=c("Necator americanus", 
                "Ancylostoma ceylanicum", "Ancylostoma duodenale")))

hw.plot=ggplot(hw, aes(x=hwepg, y=CT))+
  geom_point(aes(col=Species),alpha=0.65)+
  geom_smooth(se=FALSE,col="black",data=hw[hw$hwepg>1,], method='loess')+
  scale_y_log10(labels=yseq, breaks=yseq, limits=c(10^1.1, 10^1.7)) +
  scale_x_log10(labels=xseq, breaks=xseq, limits=c(1, 10^5))+
  xlab(expression(paste("log"[10], " Hookworm", " EPG")))+
  ylab(expression(paste("log"[10], " Hookworm Ct value")))+
  scale_color_manual(values=c(cb.blue,cb.dblue,cb.pink))+
  theme_bw()+ 
  theme(legend.position = c(0.77, 0.8), legend.background = element_rect(color = "black", 
    fill = "white", size = 0.2, linetype = "solid"))+
  ggtitle("Hookworm")+theme(plot.title = element_text(hjust = 0.5))


# tt plot
tt <- qdata.conc %>% filter(positive.Tt==1 | ttkk==1) %>%
  # impute 1 for negative values of epg
  mutate(ttepg=ifelse(ttepg==0, 1, ttepg))

tt.plot=ggplot(tt, aes(x=ttepg, y=CTmean.Tt))+
  geom_point(alpha=0.65,col=cb.green)+
  geom_smooth(se=FALSE,col="black",data=tt[tt$ttepg>1,], method='loess')+
  scale_y_log10(labels=yseq, breaks=yseq, limits=c(10^1.1, 10^1.7)) +
  scale_x_log10(labels=xseq, breaks=xseq, limits=c(1, 10^5))+
  xlab(expression(paste("log"[10], italic(" T. trichiura"), " EPG")))+
  ylab(expression(paste("log"[10], italic(" T. trichiura"), " Ct value")))+
  theme_bw()+ggtitle(expression(paste(italic("T. trichiura"))))+theme(plot.title = element_text(hjust = 0.5))

cont.plot=grid.arrange(al.plot,hw.plot,tt.plot,nrow=1)

pdf(file="~/Dropbox/WASH-B-STH-Add-on/TFGH/Results/wbb-qpcr-kk-scatter.pdf",
    width=15,height=4)
grid.draw(cont.plot)
dev.off()
