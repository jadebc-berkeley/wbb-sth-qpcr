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
# function to make pretty median (range)
# -------------------------------------
medrange=function(y, digits){
  y=y[!is.na(y)]
  min=quantile(y,prob=c(0))
  med=quantile(y,prob=c(0.5))
  max=quantile(y,prob=c(1))
  
  min=sprintf(paste0("%0.",digits,"f"),min)
  med=sprintf(paste0("%0.",digits,"f"),med)
  max=sprintf(paste0("%0.",digits,"f"),max)
  return(paste0(med, " (",min, ", ",max,")"))
}

# -------------------------------------
# KK EPG
# -------------------------------------
epg.summary <- qdata %>%
  mutate(alepg_pos = ifelse(alkk==0, NA, alepg),
         hwepg_pos = ifelse(hwkk==0, NA, hwepg),
         ttepg_pos = ifelse(ttkk==0, NA, ttepg)) %>%
  select(alepg_pos, hwepg_pos, ttepg_pos) 

epg.medrange <- sapply(epg.summary, function(x) medrange(y=x,digits=0))
epg.df <- data.frame(org = c("Ascaris lumbricoides","Hookworm","Trichuris trichiura"),
                     epg = epg.medrange)
epg.df$org = as.character(epg.df$org)

# -------------------------------------
# CT value
# -------------------------------------
CT.summary <- qdata %>%
  select(CTmean.Al, CTmean.Ac, CTmean.Ad, CTmean.Na,
         CTmean.Tt,CTmean.Ss) 

CT.medrange <- sapply(CT.summary, function(x) medrange(y=x,digits=0))

CT.summary <- matrix(t(CT.summary),ncol(CT.summary),1)

org=c("Ascaris lumbricoides","Ancylostoma ceylanicum",
      "Ancylostoma duodenale","Necator americanus","Trichuris trichiura",
      "Strongyloides stercoralis")
ct=data.frame(org=org,ct=CT.summary)
ct$org=as.character(ct$org)

tab=full_join(prev.kk,prev.q,by="org") 
tab=full_join(tab,epg.df,by="org")
tab=full_join(tab,ct,by="org")

# manual reorder
tab = tab[,c(1:3,6,4,5,7)]
tab=tab[c(1,2,4:6,3,7),]

# -------------------------------------
# save table
# -------------------------------------
write.csv(tab,file=paste0(tab_dir, "prev.table.csv"),row.names=FALSE)
