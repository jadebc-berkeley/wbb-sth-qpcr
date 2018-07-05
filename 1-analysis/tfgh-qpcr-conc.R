#######################################
# WASH Benefits Bangladesh STH KK qPCR validation
# convert CT value to DNA concentration
#######################################
library(dplyr)
library(tidyr)
library(pcr)

data.dir="/Users/jadederong/Dropbox/WASH-B-STH-Add-on/TFGH/Data/Final file/"

d=read.csv(paste0(data.dir,"Compiled controls by assay - 3-29-18 - Bangladesh_Al_nofailures.csv"))

# drop columns that take mean / SD, rename columns
d <- d %>%
  select(Plate.ID, Control.Conc, Ct.Value) %>%
  rename(id=Plate.ID, 
         conc=Control.Conc,
         ct=Ct.Value)

# not sure if we need to take the mean? 

d <- d %>%
  # group_by(id, conc) %>%
  # summarise(mean_ct=mean(ct)) %>%
  # create numeric concentration in micrograms
  mutate(conc.n=case_when(
      conc == "10ag/ul" ~ 10^(-12),
      conc == "10fg/ul" ~ 10^(-9),
      conc == "10pg/ul" ~ 10^(-6)
  )) %>%
  mutate(log.conc=log10(conc.n))

# plot concentration by mean CT value
plot(d$log.conc, d$ct)

# run linear regression
fit=glm(ct ~ log.conc, data=d)
summary(fit)

# r squared
# 1 - (Residual Deviance/Null Deviance)
1 - (976.54/12454.21)
  
# slope: 
slope=fit$coef[2]

# intercept: 
intercept=fit$coef[1]

# y = mx + b
# ct = m*(log conc) + b
# quantity of copies in a sample = 10^((ct - b)/m)
# https://www.youtube.com/watch?v=GQOnX1-SUrI

d <- d %>% mutate(copies = 10^((ct - intercept)/slope))

mean(d$copies)

#---------------------------------
# replicate using pcr package
#---------------------------------
dc=as.data.frame(d %>% select(ct))
conc=d$conc.n

pcr_assess(dc,
           amount = conc,
           method = 'standard_curve')




