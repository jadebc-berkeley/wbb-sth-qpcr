##############################################
# WASH Benefits Bangladesh STH KK qPCR validation
# between-arm comparisons
# Generate plots
##############################################
rm(list=ls())
library(ggplot2)
library(grid)
library(gridExtra)

source("~/Documents/CRG/wash-benefits/bangladesh/src/wbb-sth-qpcr/0-base-plot-functions.R")
load("~/Dropbox/WASH-B-STH-Add-on/TFGH/Data/RData/tfgh_arms.RData")

res.dir="~/Dropbox/WASH-B-STH-Add-on/TFGH/Results/"

aln=psth_n_prev_j$N.al
hwn=psth_n_prev_j$N.hw
ttn=psth_n_prev_j$N.tt
sthn=psth_n_prev_j$N.sth
alprev=al_prev
hwprev=hw_prev
ttprev=tt_prev
sthprev=sth_prev
alprh1=al_rr_h1_unadj_j
alprh2=al_rr_h2_unadj_j
hwprh1=hw_rr_h1_unadj_j
hwprh2=hw_rr_h2_unadj_j
ttprh1=tt_rr_h1_unadj_j
ttprh2=tt_rr_h2_unadj_j
sthprh1=sth_rr_h1_unadj_j
sthprh2=sth_rr_h2_unadj_j

# binary primary outcome unadj 
sth.bin.plot(
  aln=psth_n_prev_j$N.al,
  hwn=psth_n_prev_j$N.hw,
  ttn=psth_n_prev_j$N.tt,
  sthn=psth_n_prev_j$N.sth,
  alprev=al_prev,
  hwprev=hw_prev,
  ttprev=tt_prev,
  sthprev=sth_prev,
  alprh1=al_rr_h1_unadj_j,
  alprh2=al_rr_h2_unadj_j,
  hwprh1=hw_rr_h1_unadj_j,
  hwprh2=hw_rr_h2_unadj_j,
  ttprh1=tt_rr_h1_unadj_j,
  ttprh2=tt_rr_h2_unadj_j,
  sthprh1=sth_rr_h1_unadj_j,
  sthprh2=sth_rr_h2_unadj_j,
             lab="primary",res.dir=res.dir)


# epg primary outcome unadj
sth.epg.plot(psth_n_int_j$N.int.al,psth_n_int_j$N.int.hw,psth_n_int_j$N.int.tt,
             al_int_gmn,hw_int_gmn,tt_int_gmn,
             al_fecr_geo_h1_unadj_j,al_fecr_geo_h2_unadj_j,al_fecr_geo_h3_unadj_j,
             hw_fecr_geo_h1_unadj_j,hw_fecr_geo_h2_unadj_j,hw_fecr_geo_h3_unadj_j,
             tt_fecr_geo_h1_unadj_j,tt_fecr_geo_h2_unadj_j,tt_fecr_geo_h3_unadj_j,
             lab="primary")