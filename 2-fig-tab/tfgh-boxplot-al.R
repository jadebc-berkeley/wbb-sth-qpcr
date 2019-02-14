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
load("~/Dropbox/WASH-B-STH-Add-on/TFGH/Data/RData/qdata.RData")

qdata = qdata %>%
  mutate(almh.f = ifelse(almh==1, "Moderate-heavy intensity\ninfection",
                         "Light intensity\ninfection"))

pdf(file="~/Dropbox/WASH-B-STH-Add-on/TFGH/Results/wbb-alconc-boxplot.pdf",
    width=5,height=5)
ggplot(qdata,aes(y=CTmean.Al,x=almh.f))+
  geom_boxplot()+  
  geom_dotplot(aes(fill=almh.f,col=almh.f),binaxis='y',stackdir='center',
               stackratio=1,dotsize=0.7,binwidth=.008, alpha=0.5)+
  xlab("Infection intensity")+
  scale_y_log10(labels=seq(15,40,5), breaks=seq(15,40,5), limits=c(15, 40)) +
  scale_color_manual(values=c("#E69F00","#b35900"),guide=FALSE)+
  scale_fill_manual(values=c("#E69F00","#b35900"),guide=FALSE)+
  ylab(expression(paste(italic(" A. lumbricoides"), " Ct value")))+
  theme_bw()
dev.off()