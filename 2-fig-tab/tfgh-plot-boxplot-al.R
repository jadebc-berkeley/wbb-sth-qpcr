#######################################
# WASH Benefits Bangladesh STH KK qPCR validation
# boxplot of ascaris ct distribution
# within intensity categories
#######################################
rm(list=ls())
library(dplyr)
library(tidyr)
library(reshape2)
library(ggplot2)
library(grid)
library(gridExtra)
library(washb)
load("~/Dropbox/WASH-B-STH-Add-on/TFGH/Data/RData/qdata.RData")

qdata = qdata %>%
  mutate(almh.f = ifelse(almh==1, "Moderate-heavy intensity\ninfection",
                         "Light intensity\ninfection"))

pdf(file="~/Dropbox/WASH-B-STH-Add-on/TFGH/Results/wbb-alconc-boxplot.pdf",
    width=5,height=5)
ggplot(qdata,aes(y=CTmean.Al2,x=almh.f))+
  geom_boxplot()+  
  geom_dotplot(aes(fill=almh.f,col=almh.f),binaxis='y',stackdir='center',
               stackratio=1,dotsize=0.9,binwidth=.5, alpha=0.5)+
  xlab("Infection intensity")+
  scale_color_manual(values=c("#F7CA65","#BD8302"),guide=FALSE)+
  scale_fill_manual(values=c("#F7CA65","#BD8302"),guide=FALSE)+
  ylab(expression(paste(italic(" Ascaris lumbricoides"), " Cq value")))+
  theme_bw()
dev.off()


# p-value for difference in means t test
qdata = qdata %>% filter(!is.na(CTmean.Al2))

washb_ttest(Y=qdata$CTmean.Al2,tr=qdata$almh.f,strat=qdata$block, 
            contrast=c("Light intensity\ninfection","Moderate-heavy intensity\ninfection"))

# quantiles of the Ct values within infection
# intensity levels
quantile(qdata$CTmean.Al2[qdata$almh.f=="Light intensity\ninfection"], prob=c(0,.5,1))
quantile(qdata$CTmean.Al2[qdata$almh.f=="Moderate-heavy intensity\ninfection"], prob=c(0,.5,1))
