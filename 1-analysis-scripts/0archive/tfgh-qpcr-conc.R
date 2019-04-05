#######################################
# WASH Benefits Bangladesh STH KK qPCR validation
# convert CT value to DNA concentration
#######################################
library(dplyr)
library(tidyr)
library(ggplot2)

data.dir="/Users/jadederong/Dropbox/WASH-B-STH-Add-on/TFGH/Data/Final file/Compiled controls by assay - 3-29-18 - Bangladesh_"
load("~/Dropbox/WASH-B-STH-Add-on/TFGH/Data/RData/qdata.RData")

# read in datasets
al=read.csv(paste0(data.dir,"Al_nofailures.csv"))
ac=read.csv(paste0(data.dir,"Ac_nofailures.csv"))
na=read.csv(paste0(data.dir,"Na_nofailures.csv"))
ad=read.csv(paste0(data.dir,"Ad_nofailures.csv"))
tt=read.csv(paste0(data.dir,"Tt_nofailures.csv"))
ss=read.csv(paste0(data.dir,"Ss_nofailures.csv"))

al = al %>% mutate(org="al")
ac = ac %>% mutate(org="ac")
na = na %>% mutate(org="na")
ad = ad %>% mutate(org="ad")
tt = tt %>% mutate(org="tt")
ss = ss %>% mutate(org="ss")

# append datasets
d = bind_rows(al,ac,na,ad,tt,ss)

# drop columns that take mean / SD, rename columns
d <- d %>%
  select(Plate.ID, Control.Conc, Ct.Value, org) %>%
  rename(id=Plate.ID) %>%
  rename(conc=Control.Conc) %>%
  rename(ct=Ct.Value)

d <- d %>%
  # convert to attograms / ul 
  mutate(conc.n=case_when(
      conc == "10ag/ul" ~ 10^1,
      conc == "10fg/ul" ~ 10^(3),
      conc == "10pg/ul" ~ 10^(6)
  )) %>%
  mutate(log.conc=log10(conc.n))

# plot log concentration by mean CT value
ggplot(d, aes(x=log.conc,y=ct))+geom_point()+
  facet_wrap(~org)

#---------------------------------
# function for standard curve
#---------------------------------
scurve = function(log.conc, ct){
  # fit linear model
  fit = glm(ct ~ log.conc)
  
  # r squared
  r2 = 1-(fit$deviance/fit$null.deviance)
  
  # slope: 
  slope=fit$coef[2]
  
  # intercept: 
  intercept=fit$coef[1]
  
  return(data.frame(slope=slope,intercept=intercept,r2=r2))
}

org=c("al","ac","na","ad","tt","ss")

out=t(as.matrix(sapply(org,function(x) scurve(log.conc=d[d$org==x,"log.conc"],
                                                   ct=d[d$org==x,"ct"]))))
out = as.data.frame(cbind(org,out))

#---------------------------------
# function to calculate concentration
#---------------------------------
# input: intercept and slope from standard curve
# linear regression model
# y = mx + b
# ct = m*(log conc) + b
# quantity of copies in a sample = 10^((ct - b)/m)
# https://www.youtube.com/watch?v=GQOnX1-SUrI
dna.conc=function(ct, intercept, slope){
  intercept=as.numeric(intercept)
  slope=as.numeric(slope)
  return(10^((ct - intercept)/slope))
}

al <- qdata %>% 
  filter(positive.Al==1) %>%
  select(c(dataid,personid,CTmean.Al))
al <- al %>%
  mutate(copies.Al = dna.conc(ct=al$CTmean.Al, 
    intercept=out$intercept[out$org=="ac"],
    slope=out$slope[out$org=="ac"]))

ac <- qdata %>% 
  filter(positive.Ac==1) %>%
  select(c(dataid,personid,CTmean.Ac))
ac <- ac %>%
  mutate(copies.Ac = dna.conc(ct=ac$CTmean.Ac, 
    intercept=out$intercept[out$org=="ac"],
    slope=out$slope[out$org=="ac"]))

na <- qdata %>% 
  filter(positive.Na==1) %>%
  select(c(dataid,personid,CTmean.Na))
na <- na %>%
  mutate(copies.Na = dna.conc(ct=na$CTmean.Na, 
     intercept=out$intercept[out$org=="na"],
     slope=out$slope[out$org=="na"]))

ad <- qdata %>% 
  filter(positive.Ad==1) %>%
  select(c(dataid,personid,CTmean.Ad))
ad <- ad %>%
  mutate(copies.Ad = dna.conc(ct=ad$CTmean.Ad, 
     intercept=out$intercept[out$org=="ad"],
     slope=out$slope[out$org=="ad"]))

tt <- qdata %>% 
  filter(positive.Tt==1) %>%
  select(c(dataid,personid,CTmean.Tt))
tt <- tt %>%
  mutate(copies.Tt = dna.conc(ct=tt$CTmean.Tt, 
     intercept=out$intercept[out$org=="tt"],
     slope=out$slope[out$org=="tt"]))

ss <- qdata %>% 
  filter(positive.Ss==1) %>%
  select(c(dataid,personid,CTmean.Ss))
ss <- ss %>%
  mutate(copies.Ss = dna.conc(ct=ss$CTmean.Ss, 
     intercept=out$intercept[out$org=="ss"],
     slope=out$slope[out$org=="ss"]))

# merge dna concentration back on to main data
qdata.conc <- left_join(qdata, al, by=c("dataid","personid","CTmean.Al"))
qdata.conc <- left_join(qdata.conc, ac, by=c("dataid","personid","CTmean.Ac"))
qdata.conc <- left_join(qdata.conc, na, by=c("dataid","personid","CTmean.Na"))
qdata.conc <- left_join(qdata.conc, ad, by=c("dataid","personid","CTmean.Ad"))
qdata.conc <- left_join(qdata.conc, tt, by=c("dataid","personid","CTmean.Tt"))
qdata.conc <- left_join(qdata.conc, ss, by=c("dataid","personid","CTmean.Ss"))

save(qdata.conc, al, ac, na, ad, tt, ss, file="~/Dropbox/WASH-B-STH-Add-on/TFGH/Data/RData/concentration.RData")
write.csv(qdata.conc, file="~/Dropbox/WASH-B-STH-Add-on/TFGH/Data/RData/concentration.csv", row.names=FALSE)


