#######################################
# WASH Benefits Bangladesh STH KK qPCR validation
# scatter plot qpcr Cq values vs. kk EPG values
#######################################
rm(list=ls())

# configure directories, load libraries and base functions
source(paste0(here::here(), "/0-config.R"))

# load data
load(paste0(data_dir,"qdata.RData"))

#--------------------------------------
# scatter plot of qPCR CT and KK EPG results
#--------------------------------------
cb.lightorange="#E69F00"
cb.blue= "#56B4E9"
cb.green="#009E73"
cb.orange="#D55E00"
cb.pink="#CC79A7"
cb.dblue="#005787"

yseq=seq(0,45,5)
xseq=c(1, 10, 100, 1000, 10000, 100000)

# customize font size, legend
fonts = theme(legend.text.align = 0,
      axis.text.x = element_text(size=12),
      axis.title.x = element_text(size=14, margin = margin(t = 10, r = 0, b = 0, l = 0)),
      axis.text.y = element_text(size=12),
      axis.title.y = element_text(size=14, margin = margin(t = 0, r = 10, b = 0, l = 0)),
      strip.text = element_text(size=18),
      legend.text = element_text(size=10),
      legend.title = element_text(size=12),
      plot.title = element_text(hjust = 0.5))
  
# al plot
al <- qdata %>% filter(positive.Al2==1 | alkk==1) %>%
  # impute 1 for negative values of epg
  mutate(alepg=ifelse(alepg==0, 1, alepg))

al.plot=ggplot(al, aes(x=alepg, y=CTmean.Al2))+
  geom_point(alpha=0.65,col=cb.lightorange)+
  geom_smooth(se=FALSE,col="black",data=al[al$alepg>1,], method='loess')+
  scale_y_log10(labels=yseq, breaks=yseq, limits=c(10^0.5, 10^1.7)) +
  scale_x_log10(labels=xseq, breaks=xseq, limits=c(1, 10^5))+
  xlab("Kato-Katz eggs per gram")+
  ylab("qPCR Cq value")+
  theme_bw()+
  fonts + 
  ggtitle(expression(paste(italic("Ascaris lumbricoides"))))


# hw plot
hw <- qdata %>% filter(positive.Hw==1 | hwkk==1) %>%
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
  xlab("Kato-Katz eggs per gram")+
  ylab("qPCR Cq value")+
  scale_color_manual(values=c(cb.blue,cb.dblue,cb.pink))+
  theme_bw()+ 
  theme(legend.position = c(0.75, 0.8), 
        legend.background = element_rect(color = "black", 
          fill = "white", size = 0.2, linetype = "solid"),
        legend.text = element_text(face = "italic"))+
  fonts + 
  ggtitle("Hookworm")

hw.plot.poster=ggplot(hw, aes(x=hwepg, y=CT))+
  geom_point(aes(col=Species),alpha=0.65)+
  geom_smooth(se=FALSE,col="black",data=hw[hw$hwepg>1,], method='loess')+
  scale_y_log10(labels=yseq, breaks=yseq, limits=c(10^1.1, 10^1.7)) +
  scale_x_log10(labels=xseq, breaks=xseq, limits=c(1, 10^5))+
  xlab("Kato-Katz eggs per gram")+
  ylab("qPCR Ct value")+
  scale_color_manual(values=c(cb.blue,cb.dblue,cb.pink), guide=FALSE)+
  theme_bw()+ 
  ggtitle("Hookworm")+theme(plot.title = element_text(hjust = 0.5))


# tt plot
tt <- qdata %>% filter(positive.Tt==1 | ttkk==1) %>%
  # impute 1 for negative values of epg
  mutate(ttepg=ifelse(ttepg==0, 1, ttepg))

tt.plot=ggplot(tt, aes(x=ttepg, y=CTmean.Tt))+
  geom_point(alpha=0.65,col=cb.green)+
  geom_smooth(se=FALSE,col="black",data=tt[tt$ttepg>1,], method='loess')+
  scale_y_log10(labels=yseq, breaks=yseq, limits=c(10^1.2, 10^1.7)) +
  scale_x_log10(labels=xseq, breaks=xseq, limits=c(1, 10^5))+
  xlab("Kato-Katz eggs per gram")+
  ylab("qPCR Cq value")+
  theme_bw()+
  fonts + 
  ggtitle(expression(paste(italic("Trichuris trichiura"))))

cont.plot=grid.arrange(al.plot,hw.plot,tt.plot,nrow=1)

pdf(file=paste0(fig_dir,"wbb-qpcr-kk-scatter.pdf"),
    width=15,height=4)
grid.draw(cont.plot)
dev.off()

cont.plot.poster=grid.arrange(al.plot,hw.plot.poster,tt.plot,nrow=1)

pdf(file=paste0(fig_dir, "wbb-qpcr-kk-scatter-poster.pdf"),
    width=11,height=3)
grid.draw(cont.plot.poster)
dev.off()
