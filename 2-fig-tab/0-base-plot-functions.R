##############################################
# WASH Benefits Bangladesh STH Kato-Katz qPCR validation Study
# Functions to make plots
##############################################
library(ggplot2)
library(grid)
library(gridExtra)

source("~/Box Sync/WASHB Parasites/Scripts/Figures/theme_complete_bw.R")

# format point estimate and ci
pt.est.ci.f=function(obj,decimals,scale){
  a=sprintf(paste("%0.0",decimals,"f",sep=""),obj[1]*scale)
  b=sprintf(paste("%0.0",decimals,"f",sep=""),obj[2]*scale)
  c=sprintf(paste("%0.0",decimals,"f",sep=""),obj[3]*scale)
  return(paste(a," (",b,", ",c,")",sep=""))
}

# format ci
ci.f=function(obj,decimals,scale){
  b=sprintf(paste("%0.0",decimals,"f",sep=""),obj[2]*scale)
  c=sprintf(paste("%0.0",decimals,"f",sep=""),obj[3]*scale)
  return(paste("(",b,", ",c,")",sep=""))
}

# prepare prevalence data for plotting
sth.plot.prep=function(prev){
  
  prev=as.data.frame(prev)
  prev$tr=as.factor(rownames(prev))
  prev$prev.f=as.numeric(sprintf("%0.1f",prev$Prev*100))
  prev$lb.f=as.numeric(sprintf("%0.0f",prev$lb*100))
  prev$ub.f=as.numeric(sprintf("%0.0f",prev$ub*100))
  prev$x=as.factor(seq(1,4,1))
  
  return(prev)
}

# prepare prevalence data for plotting
sth.epg.plot.prep=function(mn){
  
  mn=as.data.frame(mn)
  mn$tr=as.factor(rownames(mn))
  colnames(mn)[1]="geomean"
  mn$mn.f=as.numeric(sprintf("%0.2f",mn$geomean))
  mn$x=as.factor(seq(1,7,1))
  
  return(mn)
}

# prepare prevalence data for plotting
sth.epg.em.plot.prep=function(mn){
  
  mn=as.data.frame(mn)
  mn$tr=as.factor(rownames(mn))
  mn$mn.f=as.numeric(sprintf("%0.2f",mn$geomean))
  mn$x=as.factor(seq(1,7,1))
  
  return(mn)
}









