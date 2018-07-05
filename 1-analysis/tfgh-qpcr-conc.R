#######################################
# WASH Benefits Bangladesh STH KK qPCR validation
# convert CT value to DNA concentration
#######################################
library(dplyr)
library(tidyr)
library(pcr)
library(ggplot2)

data.dir="/Users/jadederong/Dropbox/WASH-B-STH-Add-on/TFGH/Data/Final file/Compiled controls by assay - 3-29-18 - Bangladesh_"

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
  rename(id=Plate.ID, 
         conc=Control.Conc,
         ct=Ct.Value)

d <- d %>%
  mutate(conc.n=case_when(
      conc == "10ag/ul" ~ 10^(-12),
      conc == "10fg/ul" ~ 10^(-9),
      conc == "10pg/ul" ~ 10^(-6)
  )) %>%
  mutate(log.conc=log10(conc.n))

# plot concentration by mean CT value
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

al <- d %>% filter(org=="al") 
al <- al %>%
  mutate(copies = dna.conc(ct=al$ct, 
    intercept=out$intercept[out$org=="ac"],
    slope=out$slope[out$org=="ac"]))

ac <- d %>% filter(org=="ac") 
ac <- ac %>%
  mutate(copies = dna.conc(ct=ac$ct, 
    intercept=out$intercept[out$org=="ac"],
    slope=out$slope[out$org=="ac"]))

na <- d %>% filter(org=="na") 
na <- na %>%
  mutate(copies = dna.conc(ct=na$ct, 
     intercept=out$intercept[out$org=="na"],
     slope=out$slope[out$org=="na"]))

ad <- d %>% filter(org=="ad") 
ad <- ad %>%
  mutate(copies = dna.conc(ct=ad$ct, 
     intercept=out$intercept[out$org=="ad"],
     slope=out$slope[out$org=="ad"]))

tt <- d %>% filter(org=="tt") 
tt <- tt %>%
  mutate(copies = dna.conc(ct=tt$ct, 
     intercept=out$intercept[out$org=="tt"],
     slope=out$slope[out$org=="tt"]))

ss <- d %>% filter(org=="ss") 
ss <- ss %>%
  mutate(copies = dna.conc(ct=ss$ct, 
     intercept=out$intercept[out$org=="ss"],
     slope=out$slope[out$org=="ss"]))

d = bind_rows(al,ac,na,ad,tt,ss)

save(d, file="~/Dropbox/WASH-B-STH-Add-on/TFGH/Data/RData/concentration.RData")


