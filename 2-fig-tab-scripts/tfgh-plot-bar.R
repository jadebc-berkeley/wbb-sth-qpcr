#######################################
# WASH Benefits Bangladesh STH KK qPCR validation
# plot qpcr vs. kk prevalence

# by Jade Benjamin-Chung
# jadebc@berkeley.edu
#######################################
rm(list=ls())

# configure directories, load libraries and base functions
source(paste0(here::here(), "/0-config.R"))

# load data
load(paste0(data_dir,"qdata.RData"))

# drop the observation without qPCR
qdata = qdata %>% filter(!is.na(positive.Al))
assert_that(nrow(qdata)==2799)

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
  summarise_all(list(~mean(.,na.rm=TRUE))) %>%
  gather(lab,per.pos)


# calculate percent positive
n.pos = qdata %>%
  # subset columns
  select(alkk,hwkk,ttkk,
         positive.Al,positive.Tt,positive.Na,positive.Ad,
         positive.Ac) %>%
  # calculate sum of positives
  summarise_all(list(~sum(.,na.rm=TRUE)))  %>%
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


pdf(file=paste0(fig_dir, "wbb-qpcr-kk-bargraph.pdf"),
    width=10,height=4)
ggplot(bar.data,aes(x=test,y=per.pos,fill=org),col="black")+
  geom_bar(aes(fill=org),stat="identity",colour="black",
           position='stack') + 
  facet_grid(~ orgcat)+
  scale_fill_manual("Organism",values=mycol,
                    labels=expression(italic("Ascaris lumbricoides"),
                                      "Hookworm",
                                      italic("Ancylostoma ceylanicum"),
                                      italic("Ancylostoma duodenale"),
                                      italic("Necator americanus"),
                                      italic("Trichuris trichiura"))) +
  scale_y_continuous(limits=c(0,40)) +
  geom_text(aes(label=per.f,
                vjust=c(-0.3,-0.3,-0.3,-0.3,-0.3,-2.2,-0.3,-0.3)), size=5) +
  theme_bw() +
  
  # customize font size, legend
  theme(legend.text.align = 0,
        axis.text.x = element_text(size=12),
        axis.title.x = element_text(size=14, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size=14, margin = margin(t = 0, r = 10, b = 0, l = 0)),
        strip.text = element_text(size=14),
        legend.text = element_text(size=12),
        legend.title = element_text(size=14)) +
  
  ylab("Prevalence") +
  xlab("Diagnostic method")
dev.off()

pdf(file=paste0(fig_dir, "wbb-qpcr-kk-bargraph-poster.pdf"),
    width=4,height=3.5)
ggplot(bar.data,aes(x=test,y=per.pos,fill=org),col="black")+
  geom_bar(aes(fill=org),stat="identity",colour="black",
           position='stack')+facet_grid(~orgcat)+
  scale_fill_manual("Organism",values=mycol) +
  scale_y_continuous(limits=c(0,40))+
  geom_text(aes(label=per.f,vjust=c(-0.3,-0.3,-0.3,-0.3,-0.3,-1.8,-0.3,-0.3)))+
  theme_bw()+ylab("Prevalence")+xlab("Diagnostic method")+
  theme(legend.position="bottom")
dev.off()

