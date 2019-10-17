#######################################
# WASH Benefits Bangladesh STH KK qPCR validation
# Clean and merge qPCR and KK data

# by Jade Benjamin-Chung
# jadebc@berkeley.edu
#######################################
rm(list=ls())

# configure directories, load libraries and base functions
source(paste0(here::here(), "/0-config.R"))

#--------------------------------------
# read in qPCR data
#--------------------------------------
raw.data.dir="~/Dropbox/WASH-B-STH-Add-on/TFGH/Data/Final file/"
res.name="FINAL FINAL FINAL Master Results File - Bangladesh - 4-30-18"
iac=read.csv(file=paste0(raw.data.dir,res.name,"-IAC.csv"),stringsAsFactors=FALSE)
na=read.csv(file=paste0(raw.data.dir,res.name,"-NA.csv"),stringsAsFactors=FALSE)
ad=read.csv(file=paste0(raw.data.dir,res.name,"-AD.csv"),stringsAsFactors=FALSE)
ss=read.csv(file=paste0(raw.data.dir,res.name,"-SS.csv"),stringsAsFactors=FALSE)
tt=read.csv(file=paste0(raw.data.dir,res.name,"-TT.csv"),stringsAsFactors=FALSE)
al=read.csv(file=paste0(raw.data.dir,res.name,"-AL.csv"),stringsAsFactors=FALSE)
ac=read.csv(file=paste0(raw.data.dir,res.name,"-AC.csv"),stringsAsFactors=FALSE)

colnames=c("sampleid","sampleno","assay",
  "CTmean","CTSD","assaydate","plateid","comments")

colnames(iac)=colnames
colnames(na)=colnames
colnames(ad)=colnames
colnames(ss)=colnames
colnames(tt)=colnames
colnames(al)=colnames
colnames(ac)=colnames

# drop problematic sample
na=na[na$comments=="",]
ad=ad[ad$comments=="",]
ss=ss[ss$comments=="",]
tt=tt[tt$comments=="",]
ac=ac[ac$comments=="",]

# correcting assay column in Na
na <- na %>% mutate(assay="Na")

qpcr=rbind(iac,na,ad,ss,tt,al,ac)
qpcr=qpcr[!is.na(qpcr$sampleid),]

# check number of results per assay
table(qpcr$assay)

# correcting data entry error
qpcr$sampleid[qpcr$sampleid==":559901ETS"]="559901ETS1"

# take mean of results for each sample
qpcr.mean <- qpcr %>%
  group_by(sampleid,assay) %>%
  summarise(CTmean=mean(CTmean),CTSD=mean(CTSD)) %>%
  mutate(positive=ifelse(CTmean<40,1,0)) %>%
  mutate(positive=ifelse(is.na(CTmean),0,positive))

mean.l <- qpcr.mean %>% 
  group_by(sampleid) %>%
  dplyr::select(sampleid,assay,CTmean) %>% 
  spread(assay,CTmean) %>%
  rename(CTmean.Ac=Ac, CTmean.Ad=Ad, CTmean.Al=Al,
         CTmean.IAC=IAC, CTmean.Na=Na, 
         CTmean.Ss=Ss,CTmean.Tt=Tt)

sd.l <- qpcr.mean %>% 
  group_by(sampleid) %>%
  dplyr::select(sampleid,assay,CTSD) %>% 
  spread(assay,CTSD) %>%
  rename(CTSD.Ac=Ac, CTSD.Ad=Ad, CTSD.Al=Al,
         CTSD.IAC=IAC, CTSD.Na=Na, 
         CTSD.Ss=Ss,CTSD.Tt=Tt)

pos.l <- qpcr.mean %>% 
  group_by(sampleid) %>%
  dplyr::select(sampleid,assay,positive) %>% 
  spread(assay,positive) %>%
  rename(positive.Ac=Ac, positive.Ad=Ad, positive.Al=Al,
         positive.IAC=IAC, positive.Na=Na, 
         positive.Ss=Ss,positive.Tt=Tt) 

qpcr.w <- full_join(mean.l, sd.l, by=c("sampleid"))
qpcr.w <- full_join(qpcr.w, pos.l, by=c("sampleid"))

qpcr.w$dataid=substr(qpcr.w$sampleid,1,5)
qpcr.w$personid=paste(substr(qpcr.w$sampleid,7,7),1,sep="")
qpcr.w$qpcr="Done"

qpcr.w <- qpcr.w %>%
  mutate(positive.Hw=case_when(
    positive.Ac==1 | positive.Na==1 | positive.Ad==1 ~ 1,
    positive.Ac==0 & positive.Na==0 & positive.Ad==0 ~ 0,
    TRUE ~ NA_real_
  )) %>%
  # manual correction of id that doesn't match list sent to Smith
  ungroup() %>%
  mutate(sampleid=ifelse(sampleid=="559901ETS1","59901ETS1",sampleid)) %>%
  mutate(sampleid=ifelse(sampleid=="18705EOS1","18705ECS1",sampleid)) %>%
  mutate(personid=substr(sampleid,7,7),
         dataid=substr(sampleid,1,5)) %>%
  mutate(personid=ifelse(personid=="T","T1",personid)) %>%
  mutate(personid=ifelse(personid=="C","C1",personid)) %>%
  mutate(personid=ifelse(personid=="O","O1",personid)) %>%
  mutate(personid=ifelse(personid=="W","T2",personid)) 

# confirm that sample 10805ETS1 was dropped from qPCR 
qpcr.w = qpcr.w %>% 
  mutate(positive.Al = ifelse(sampleid=="10805ETS1", NA, positive.Al))

assert_that(all(is.na(qpcr.w %>% 
              filter(sampleid=="10805ETS1") %>%
                dplyr::select(CTmean.Ac, CTmean.Ad, CTmean.Al, CTmean.Na, CTmean.Ss, CTmean.Tt,
                     CTSD.Ac, CTSD.Ad, CTSD.Al, CTSD.Na, CTSD.Ss, CTSD.Tt,
                     positive.Ac, positive.Ad, positive.Al, positive.Na, positive.Ss, positive.Tt))))

#--------------------------------------
# read in revised Ascaris qPCR data
#--------------------------------------
ascaris_new=read.csv("~/Dropbox/WASH-B-STH-Add-on/TFGH/Data/Revised ascaris assay/KK v New assay comparison_4-3-19_Updated with data for missing samples.csv", stringsAsFactors=FALSE)
nrow(ascaris_new)

colnames(ascaris_new) = c("sampleid", "alepg", "al", "X", "CTmean.Al2", "CTSD.Al2", "XX")
ascaris_new$dataid=substr(ascaris_new$sampleid,1,5)
ascaris_new$personid=paste(substr(ascaris_new$sampleid,7,7),1,sep="")
ascaris_new = ascaris_new %>% dplyr::select(-c(X,XX, alepg, al))

# indicator for positive Al 
ascaris_new <- ascaris_new %>%
  mutate(CTmean.Al2 = as.numeric(CTmean.Al2)) %>%
  mutate(positive.Al2=ifelse(CTmean.Al2<40,1,0)) %>%
  mutate(positive.Al2=ifelse(is.na(CTmean.Al2),0,positive.Al2))

# manual correction of id that doesn't match list sent to Smith
ascaris_new = ascaris_new %>%
  ungroup() %>%
    mutate(sampleid=ifelse(sampleid=="18705EOS1","18705ECS1",sampleid))

# manual correction of data for id missing in xls file
# Nils sent based on Nils' email 
ascaris_new %>% filter(sampleid=="59901ETS1")

addid = data.frame(
  sampleid="59901ETS1",
  CTmean.Al2 = NA,
  CTSD.Al2 = NA,
  dataid = "59901",
  personid = "T1",
  positive.Al2 = 0
)

ascaris_new = bind_rows(ascaris_new, addid)

# create standard ids
ascaris_new = ascaris_new %>%
    mutate(personid=substr(sampleid,7,7),
           dataid=substr(sampleid,1,5)) %>%
    mutate(personid=ifelse(personid=="T","T1",personid)) %>%
    mutate(personid=ifelse(personid=="C","C1",personid)) %>%
    mutate(personid=ifelse(personid=="O","O1",personid)) %>%
    mutate(personid=ifelse(personid=="W","T2",personid)) 


# confirm that sample 10805ETS1 was dropped from qPCR 
ascaris_new = ascaris_new %>% 
  mutate(positive.Al2 = ifelse(sampleid=="10805ETS1", NA, positive.Al2))

assert_that(all(is.na(ascaris_new %>% 
                        filter(sampleid=="10805ETS1") %>%
                        dplyr::select(CTmean.Al2, CTSD.Al2, positive.Al2))))

# # confirm that only one sample has a missing ascaris result
assert_that(length(ascaris_new$positive.Al2[is.na(ascaris_new$positive.Al2)]) == 1)

#--------------------------------------
# merge in kk data
#--------------------------------------
kk=read.csv("~/Dropbox/WASHB Parasites/Analysis datasets/Jade/sth.csv")
kk$dataid=as.character(kk$dataid)
kk$personid=as.character(kk$personid)
kk=kk[,c("dataid","personid","block","clusterid","tr",
         "sex","dw","aged","agem","agey","counter1","counter2","labdate",
      "alepg","hwepg","ttepg",
       "logalepg","loghwepg","logttepg",
       "al","tt","hw","sth","alint","ttint","hwint")]

colnames(kk)[which(colnames(kk)=="tt")]="ttkk"
colnames(kk)[which(colnames(kk)=="al")]="alkk"
colnames(kk)[which(colnames(kk)=="hw")]="hwkk"

# merge kk and qPCR data
data=full_join(kk,qpcr.w,by=c("dataid","personid"))
data$qpcr[is.na(data$qpcr)]="Not done"

# subset to rows in which qpcr was done
qdata=data[data$qpcr=="Done",]
nrow(qdata)

# merge on re-run ascaris qPCR data
qdata=full_join(qdata,ascaris_new,by=c("dataid","personid","sampleid"))
nrow(qdata)

# confirm that all rows match between 
# first and second batch of assays
anti = anti_join(qdata,ascaris_new,by=c("dataid","personid","sampleid"))
assert_that(nrow(anti)==0)

# -------------------------------------------
# create gold standard by pooling together 
# kk and qPCR results
# -------------------------------------------
qdata = qdata %>% 
  mutate(positive.Hw=ifelse(positive.Na==1 | positive.Ac==1 | positive.Ad==1,1,0)) %>%
  mutate(gold.hwpos=ifelse(hwkk==1 | positive.Hw==1,1,0)) %>%
  mutate(gold.ttpos=ifelse(ttkk==1 | positive.Tt==1,1,0)) %>%
  mutate(gold.alpos=ifelse(alkk==1 | positive.Al2==1,1,0)) %>%
  mutate(gold.sthpos=ifelse(sth==1 | positive.Al2==1| positive.Hw==1 | positive.Tt==1,1,0))

# -------------------------------------------
# create indicators for moderate/heavy kk infection
# -------------------------------------------
qdata <- qdata %>%
  mutate(
    almh = case_when(alepg >= 5000 & alkk == 1 ~ 1,
                     alepg < 5000 & alkk == 1 ~ 0,
                     alkk == 0 ~ NA_real_),
    hwmh = case_when(hwepg >= 2000 & hwkk == 1 ~ 1,
                     hwepg < 2000 & hwkk == 1 ~ 0,
                     hwkk == 0 ~ NA_real_),
    ttmh = case_when(ttepg >= 1000 & ttkk == 1 ~ 1,
                     ttepg < 1000 & ttkk == 1 ~ 0,
                     ttkk == 0 ~ NA_real_),
    
    al_light = case_when(alepg < 5000 & alkk == 1 ~ 1,
                         alepg >= 5000 & alkk == 1 ~ 0,
                         alkk == 0 ~ NA_real_),
    hw_light = case_when(hwepg < 2000 & hwkk == 1 ~ 1,
                         hwepg >= 2000 & hwkk == 1 ~ 0,
                         hwkk == 0 ~ NA_real_),
    tt_light = case_when(ttepg < 1000 & ttkk == 1 ~ 1,
                         ttepg >= 1000 & ttkk == 1 ~ 0,
                         ttkk == 0 ~ NA_real_),
    
    al_mod = case_when( alepg >= 5000 & alepg < 50000 & alkk == 1 ~ 1,
                       (alepg < 5000 | alepg >= 50000) & alkk == 1 ~ 0,
                        alkk == 0 ~ NA_real_),
    hw_mod = case_when( hwepg >= 2000 & hwepg < 4000 & hwkk == 1 ~ 1,
                       (hwepg < 2000 | hwepg >= 4000) & hwkk == 1 ~ 0,
                        hwkk == 0 ~ NA_real_),
    tt_mod = case_when( ttepg >= 1000 & ttepg < 10000 & ttkk == 1 ~ 1,
                       (ttepg < 1000 | ttepg >= 10000) & ttkk == 1 ~ 0,
                        ttkk == 0 ~ NA_real_),
    
    al_heavy = case_when(alepg >= 50000 & alkk == 1 ~ 1,
                         alepg < 50000 & alkk == 1 ~ 0,
                         alkk == 0 ~ NA_real_),
    hw_heavy = case_when(hwepg >= 4000 & hwkk == 1 ~ 1,
                         hwepg < 4000 & hwkk == 1 ~ 0,
                         hwkk == 0 ~ NA_real_),
    tt_heavy = case_when(ttepg >= 10000 & ttkk == 1 ~ 1,
                         ttepg < 10000 & ttkk == 1 ~ 0,
                         ttkk == 0 ~ NA_real_)) %>%
  mutate(
    almh.f = case_when(almh == 1 ~ "Moderate-heavy intensity",
                       almh == 0 ~ "Light intensity",
                       is.na(almh) ~ "Uninfected"),
    
    hwmh.f = case_when(hwmh == 1 ~ "Moderate-heavy intensity",
                       hwmh == 0 ~ "Light intensity",
                       is.na(hwmh) ~ "Uninfected"),
    
    ttmh.f = case_when(ttmh == 1 ~ "Moderate-heavy intensity",
                       ttmh == 0 ~ "Light intensity",
                       is.na(ttmh) ~ "Uninfected"),
    
    al_light.f = case_when(al_light == 1 ~ "Light intensity",
                           al_light == 0 ~ "Moderate-heavy intensity",
                           is.na(al_light) ~ "Uninfected"),
    
    hw_light.f = case_when(hw_light == 1 ~ "Light intensity",
                           hw_light == 0 ~ "Moderate-heavy intensity",
                           is.na(hw_light) ~ "Uninfected"),
    
    tt_light.f = case_when(tt_light == 1 ~ "Light intensity",
                           tt_light == 0 ~ "Moderate-heavy intensity",
                           is.na(tt_light) ~ "Uninfected"),
    
    al_mod.f = case_when(al_mod == 1 ~ "Moderate intensity",
                         al_mod == 0 ~ "Light or heavy intensity",
                         is.na(al_mod) ~ "Uninfected"),
    
    hw_mod.f = case_when(hw_mod == 1 ~ "Moderate intensity",
                         hw_mod == 0 ~ "Light or heavy intensity",
                         is.na(hw_mod) ~ "Uninfected"),
    
    tt_mod.f = case_when(tt_mod == 1 ~ "Moderate intensity",
                         tt_mod == 0 ~ "Light or heavy intensity",
                         is.na(tt_mod) ~ "Uninfected"),
    
    al_heavy.f = case_when(al_heavy == 1 ~ "Heavy intensity",
                           al_heavy == 0 ~ "Light-moderate intensity",
                           is.na(al_heavy) ~ "Uninfected"),
    
    hw_heavy.f = case_when(hw_heavy == 1 ~ "Heavy intensity",
                           hw_heavy == 0 ~ "Light-moderate intensity",
                           is.na(hw_heavy) ~ "Uninfected"),
    
    tt_heavy.f = case_when(tt_heavy == 1 ~ "Heavy intensity",
                           tt_heavy == 0 ~ "Light-moderate intensity",
                           is.na(tt_heavy) ~ "Uninfected")
  )

# -------------------------------------------
# create any STH
# -------------------------------------------
qdata = qdata %>%
  mutate(positive.Sth = ifelse(
      positive.Al2 == 1 |
      positive.Na == 1 |
      positive.Ac == 1 | 
      positive.Ad == 1 |
      positive.Tt == 1 , 1,  0
  ))

# -------------------------------------------
# drop original Al assay, replace with new assay
# -------------------------------------------
qdata = qdata %>% 
  dplyr::select(-c(positive.Al, CTmean.Al, CTSD.Al)) %>%
  rename(positive.Al = positive.Al2,
         CTmean.Al = CTmean.Al2,
         CTSD.Al = CTSD.Al2) 

qdata = qdata %>%
  dplyr::select(clusterid, dataid, block,  personid, sampleid, tr, sex, dw,
         aged,   agem,   agey,   
         counter1, counter2,
         alepg,  hwepg,  ttepg, 
         logalepg, loghwepg, logttepg, alkk,   ttkk,   hwkk,  
         sth, alint,  ttint,  hwint,  
         
         almh, hwmh, ttmh, 
         almh.f, hwmh.f, ttmh.f,
         
         al_light, hw_light, tt_light,
         al_light.f, hw_light.f, tt_light.f,
         
         al_mod, hw_mod, tt_mod,
         al_mod.f, hw_mod.f, tt_mod.f,
         
         al_heavy, hw_heavy, tt_heavy,
         al_heavy.f, hw_heavy.f, tt_heavy.f, 
         
         CTmean.Ac, CTmean.Ad,
         CTmean.Al, CTmean.IAC,   CTmean.Na, CTmean.Ss, CTmean.Tt, CTSD.Ac, 
         CTSD.Ad,  CTSD.Al,  CTSD.IAC, CTSD.Na,  CTSD.Ss,  CTSD.Tt, 
         positive.Ac,  positive.Ad,  positive.Al,  positive.IAC, positive.Na,  positive.Ss,
         positive.Tt, positive.Hw, positive.Sth, everything())

# confirm n = 2800
assert_that(nrow(qdata)==2800)

#--------------------------------------
# save data
#--------------------------------------
# qpcr data
save(qdata,file=paste0(data_dir ,"qdata.RData"))
write.csv(qdata,file=paste0(data_dir, "qdata.csv"),row.names=FALSE)

# all kk and partial qpcr data 
save(data,file=paste0(data_dir ,"data.RData"))
saveRDS(data,file=paste0(data_dir ,"data.RDS"))
write.csv(data,file=paste0(data_dir, "data.csv"),row.names=FALSE)


