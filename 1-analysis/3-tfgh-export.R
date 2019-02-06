# export kk data for steve w

load("~/Dropbox/WASH-B-STH-Add-on/TFGH/Data/RData/qdata.RData")
library(dplyr)


out = qdata %>%
  select(dataid, personid, alepg, hwepg, ttepg, alkk, hwkk, ttkk) %>%
  mutate(sampleid=paste0(dataid,"E",substr(personid,1,1),
                         "S1")) %>%
  select(-c(dataid,personid)) %>%
  select(sampleid,everything())

nrow(out)

write.csv(out, file="~/Dropbox/WASH-B-STH-Add-on/TFGH/Data/wbb-kk.csv",
          row.names=FALSE)