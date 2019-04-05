#######################################
# WASH Benefits Bangladesh STH KK qPCR validation
# boxplot of ascaris ct distribution
# within intensity categories

# by Jade Benjamin-Chung
# jadebc@berkeley.edu
#######################################
rm(list=ls())

# configure directories, load libraries and base functions
source(paste0(here::here(), "/0-config.R"))

load("~/Dropbox/WASH-B-STH-Add-on/TFGH/Data/RData/qdata.RData")

qdata = qdata %>%
  mutate(almh.f = ifelse(almh==1, "Moderate-heavy intensity\ninfection",
                         "Light intensity\ninfection"))

pdf(file=paste0(fig_dir,"wbb-alconc-boxplot.pdf"),
    width=5,height=5)
ggplot(qdata,aes(y=CTmean.Al2,x=almh.f))+
  geom_boxplot()+  
  geom_dotplot(aes(fill=almh.f,col=almh.f),binaxis='y',stackdir='center',
               stackratio=1,dotsize=0.9,binwidth=.5, alpha=0.5)+
  xlab("Infection intensity")+
  scale_color_manual(values=c("#F7CA65","#BD8302"),guide=FALSE)+
  scale_fill_manual(values=c("#F7CA65","#BD8302"),guide=FALSE)+
  ylab(expression(paste(italic(" Ascaris lumbricoides"), " Cq value")))+
  theme_bw() +

# customize font size, legend
theme(legend.text.align = 0,
      axis.text.x = element_text(size=12),
      axis.title.x = element_text(size=14, margin = margin(t = 10, r = 0, b = 0, l = 0)),
      axis.text.y = element_text(size=12),
      axis.title.y = element_text(size=14, margin = margin(t = 0, r = 10, b = 0, l = 0))) 
  
dev.off()


# p-value for difference in means t test
qdata = qdata %>% filter(!is.na(CTmean.Al2))

washb_ttest(Y=qdata$CTmean.Al2,tr=qdata$almh.f,strat=qdata$block, 
            contrast=c("Light intensity\ninfection","Moderate-heavy intensity\ninfection"))

# quantiles of the Ct values within infection
# intensity levels
quantile(qdata$CTmean.Al2[qdata$almh.f=="Light intensity\ninfection"], prob=c(0,.5,1))
quantile(qdata$CTmean.Al2[qdata$almh.f=="Moderate-heavy intensity\ninfection"], prob=c(0,.5,1))
