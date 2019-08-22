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

n_modheav = length(qdata$almh[qdata$almh==1 & qdata$alkk==1])
n_light = length(qdata$almh[qdata$almh==0 & qdata$alkk==1])
n_negkk = length(qdata$alkk[qdata$alkk==0])

assert_that(n_modheav + n_light + n_negkk <=2800) 

qdata = qdata %>%
  mutate(xcat = case_when(
    almh==1 & alkk==1 ~ paste0("Moderate-heavy\nintensity\ninfection\n\n(N=",n_modheav,")"),
    almh==0 & alkk==1 ~ paste0("Light intensity\ninfection\n\n\n(N=", n_light, ")"),
    alkk==0 ~ paste0("Kato-Katz\nNegative\n\n\n(N=",n_negkk,")")
    )
    )

darkgreen = "#026813"
green = "#5DE674"
red = "#FA2A27"

pdf(file=paste0(fig_dir,"wbb-alconc-boxplot.pdf"),
    width=5,height=5)
ggplot(qdata,aes(y=CTmean.Al,x=xcat))+
  geom_dotplot(aes(fill=xcat,col=xcat),binaxis='y',stackdir='center',
               stackratio=1,dotsize=0.9,binwidth=.5, alpha=0.5)+
  geom_boxplot(alpha = 0)+  
  xlab("")+
  scale_color_manual(values=c(red, green, darkgreen),guide=FALSE)+
  scale_fill_manual(values=c(red, green, darkgreen),guide=FALSE)+
  ylab(expression(paste(italic(" Ascaris lumbricoides"), " Cq value")))+
  theme_bw() +
  
# customize font size, legend
  theme(
      axis.text.x = element_text(size=12),
      axis.title.x = element_text(size=14, margin = margin(t = 10, r = 0, b = 0, l = 0)),
      axis.text.y = element_text(size=12),
      axis.title.y = element_text(size=14, margin = margin(t = 0, r = 10, b = 0, l = 0)),
      legend.position = "none") 
  
dev.off()


# p-value for difference in means t test
qdata = qdata %>% filter(!is.na(CTmean.Al)) %>%
  mutate(almh.f = case_when(
    almh==1 & alkk==1 ~ "Moderate-heavy intensity infection",
    almh==0 & alkk==1 ~ "Light intensity infection",
    alkk==0 ~ "Kato-Katz negative"
  ))

washb_ttest(Y=qdata$CTmean.Al,tr=qdata$almh.f,strat=qdata$block, 
            contrast=c("Light intensity infection",
                       "Moderate-heavy intensity infection"))

washb_ttest(Y=qdata$CTmean.Al,tr=qdata$almh.f,strat=qdata$block, 
            contrast=c("Light intensity infection",
                       "Kato-Katz negative"))

washb_ttest(Y=qdata$CTmean.Al,tr=qdata$almh.f,strat=qdata$block, 
            contrast=c("Moderate-heavy intensity infection",
                       "Kato-Katz negative"))


# quantiles of the Ct values within infection
# intensity levels
quantile(qdata$CTmean.Al[qdata$almh.f=="Light intensity infection"], prob=c(0,.5,1))
quantile(qdata$CTmean.Al[qdata$almh.f=="Moderate-heavy intensity infection"], prob=c(0,.5,1))
quantile(qdata$CTmean.Al[qdata$almh.f=="Kato-Katz negative"], prob=c(0,.5,1))
