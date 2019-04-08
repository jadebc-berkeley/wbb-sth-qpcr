#######################################
# WASH Benefits Bangladesh STH KK qPCR validation
# table of qpcr vs. kk prevalence and intensity

# by Jade Benjamin-Chung
# jadebc@berkeley.edu
#######################################
rm(list=ls())

# configure directories, load libraries and base functions
source(paste0(here::here(), "/0-config.R"))

# load data
load(paste0(data_dir,"qdata.RData"))
load(paste0(data_dir,"prev_results.RData"))

# -------------------------------------
# prevalence
# -------------------------------------
prev.kk=as.data.frame(rbind(al.kk,hw.kk,tt.kk))
prev.kk$org=c("Ascaris lumbricoides","Hookworm","Trichuris trichiura")
colnames(prev.kk)[5:6]=c("lb","ub")
prev.kk = prev.kk %>%
  mutate(prevkk=ptestci.format(Mean,lb,ub,decimals=1,scale=100)) %>%
  select(org,Nkk=N,prevkk) 

prev.q=as.data.frame(rbind(al.q,hw.q,na.q,ac.q,ad.q,tt.q,ss.q))
prev.q$org=c("Ascaris lumbricoides","Hookworm",
           "Necator americanus","Ancylostoma ceylanicum","Ancylostoma duodenale",
           "Trichuris trichiura","Strongyloides stercoralis")
colnames(prev.q)[5:6]=c("lb","ub")
prev.q = prev.q %>%
  mutate(prevq=ptestci.format(Mean,lb,ub,decimals=1,scale=100)) %>%
  select(org,Nq=N,prevq)

# -------------------------------------
# KK EPG
# -------------------------------------
epg=as.data.frame(rbind(al.kk.gmn,hw.kk.gmn,tt.kk.gmn))
epg$org=c("Ascaris lumbricoides","Hookworm","Trichuris trichiura")
colnames(epg)[2:3]=c("lb","ub")
epg = epg %>%
  mutate(epg=ptestci.format(Mean,lb,ub,decimals=2,scale=1)) %>%
  select(org,epg)

# -------------------------------------
# CT value
# -------------------------------------
# function to make pretty median (range)
medrange=function(y){
  y=y[!is.na(y)]
  min=quantile(y,prob=c(0))
  med=quantile(y,prob=c(0.5))
  max=quantile(y,prob=c(1))
  
  min=sprintf("%0.1f",min)
  med=sprintf("%0.1f",med)
  max=sprintf("%0.1f",max)
  return(paste0(med, " (",min, ", ",max,")"))
}

CT.summary <- qdata %>%
  select(CTmean.Al, CTmean.Ac, CTmean.Ad, CTmean.Na,
         CTmean.Tt,CTmean.Ss) %>%
  summarise_all(list(medrange)) 

CT.summary <- matrix(t(CT.summary),ncol(CT.summary),1)

org=c("Ascaris lumbricoides","Ancylostoma ceylanicum",
      "Ancylostoma duodenale","Necator americanus","Trichuris trichiura",
      "Strongyloides stercoralis")
ct=data.frame(org=org,ct=CT.summary)
ct$org=as.character(ct$org)

tab=full_join(prev.kk,prev.q,by="org") 
tab=full_join(tab,epg,by="org")
tab=full_join(tab,ct,by="org")

# manual reorder
tab = tab[,c(1:3,6,4,5,7)]
tab=tab[c(1,2,4:6,3,7),]

# -------------------------------------
# save table
# -------------------------------------
write.csv(tab,file=paste0(tab_dir, "prev.table.csv"),row.names=FALSE)
