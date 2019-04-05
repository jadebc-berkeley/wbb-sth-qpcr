# make single plot
makeplot=function(n,prev,prh1,prh2,ytitle,ylim){
  
  # format estimates
  n=paste("(N=",format(n,big.mark=","),")",sep="")
  rr=apply(as.matrix(prh1),1,pt.est.ci.f,decimals=2,scale=1)
  rrh1=c("ref",rr)
  rrh2w=c("","ref","",pt.est.ci.f(prh2[1,],decimals=2,scale=1))
  rrh2s=c("","","ref",pt.est.ci.f(prh2[2,],decimals=2,scale=1))
  
  cliney=as.numeric(prev$prev.f[1])
  
  # define color palette
  black = "#000004FF"
  blue = "#3366AA"
  teal = "#11AA99"
  green = "#66AA55"
  chartr = "#CCCC55"
  magent = "#992288"
  red = "#EE3333"
  orange = "#EEA722"
  yellow = "#FFEE33"
  grey = "#777777"
  cols=c(black,blue,teal,orange)
  
  colh2w=c("#5e5f60",cols[2],"#5e5f60","#5e5f60")
  colh2s=c("#5e5f60","#5e5f60",cols[3],"#5e5f60")
  # colh2h=c("#5e5f60","#5e5f60","#5e5f60",cols[4],"#5e5f60","#5e5f60","#5e5f60")
  # colh3wsh=c("#5e5f60","#5e5f60","#5e5f60","#5e5f60",cols[5],"#5e5f60","#5e5f60")
  # colh3n=c("#5e5f60","#5e5f60","#5e5f60","#5e5f60","#5e5f60",cols[6],"#5e5f60")
  
  if(ytitle=="Ascaris"){
    ytitlex=-0.75
    ylim=c(0,60)
    pylim=58
    linediff=0.038
    secdiff=0.07  
    gap=5
  }
  if(ytitle=="Hookworm"){
    ytitlex=-0.6
    ylim=c(0,20)
    pylim=20
    linediff=0.06
    secdiff=0.1
    gap=2
  }
  if(ytitle=="Trichuris"){
    ytitlex=-0.7
    ylim=c(0,20)
    pylim=20
    linediff=0.06
    secdiff=0.1
    gap=2
  }
  if(ytitle=="Any STH"){
    ytitlex=-0.69
    ylim=c(20,60)
    pylim=58
    linediff=0.038
    secdiff=0.07  
    gap=5
  }
  
  g1=ggplot(prev,aes(x=x,y=prev.f))+
    geom_point(aes(col=x),alpha=0.7,size=1.5,show.legend=FALSE)+
    geom_errorbar(aes(ymin=lb.f,ymax=ub.f,col=x),width=0.11,show.legend=FALSE)+
    geom_hline(yintercept=cliney,linetype="dashed")+
    coord_cartesian(ylim = ylim,xlim=c(1,4)) +
    scale_color_manual("",values=cols)+
    ylab("Prevalence\nat 2-year\nfollow-up (%)")+xlab("")+
    scale_y_continuous(breaks=seq(ylim[1],ylim[2],gap),labels=seq(ylim[1],ylim[2],gap))+
    theme_complete_bw() +
    theme(plot.margin = unit(c(6.8, 0.5, 0.5, 3.5), "lines"))+
    
    annotate(geom="text",x=0.65,y=cliney+0.05,label="Control",size=2.5)+
    annotate(geom="text",x=0.62,y=cliney+0.02,label="mean",size=2.5)+
    annotate(geom="text",x=ytitlex,y=pylim*(1.1+linediff*4+secdiff*2.95),label=ytitle,size=5.5)+
    annotate(geom="text",x=seq(1,4,1),y=pylim*(1.1+linediff*4+secdiff*3.4),
             label=c("Control","Water","Sanitation","Combined"),
             color=cols[1:4],size=3)+
    annotate(geom="text",x=seq(1,4,1),y=pylim*(1.1+linediff*4+secdiff*2.8),
             label=c("","","","WSH"),
             color=cols[1:4],size=3)+
    annotate(geom="text",x=seq(1,4,1),y=pylim*(1.1+linediff*4+secdiff*2.1),
             label=n,size=3,color="#5e5f60")+
    annotate(geom="text",x=-0.3,y=pylim*(1.1+linediff*4+secdiff*2.02), label="Prevalence Ratio (95% CI)",size=3.5)+
    annotate(geom="text",x=-0.1,y=pylim*(1.1+linediff*3+secdiff*2), label="Intervention vs. Control",size=2.7,col="#5e5f60")+
    
    annotate(geom="text",x=seq(1,4,1),y=pylim*(1.1+linediff*3+secdiff*2),label=rrh1,size=2.7,col="#5e5f60")+
    
    annotate(geom="text",x=0.22,y=pylim*(1.1+linediff*3+secdiff), label="WSH vs. W",size=2.7,col="#5e5f60")+
    annotate(geom="text",x=0.22,y=pylim*(1.1+linediff*2+secdiff), label="WSH vs. S",size=2.7,col="#5e5f60")+
    # annotate(geom="text",x=0.22,y=pylim*(1.1+linediff+secdiff), label="WSH vs. H",size=2.7,col="#5e5f60")+
    annotate(geom="text",x=seq(1,4,1),y=pylim*(1.1+linediff*3+secdiff),label=rrh2w,size=2.7,col=colh2w)+
    annotate(geom="text",x=seq(1,4,1),y=pylim*(1.1+linediff*2+secdiff),label=rrh2s,size=2.7,col=colh2s)+
    # annotate(geom="text",x=seq(1,4,1),y=pylim*(1.1+linediff+secdiff),label=rrh2h,size=2.7,col=colh2h)+
    
    annotate(geom="text",x=seq(1,4,1)+.15,y=prev$prev.f,label=prev$prev.f,size=2.7,col=cols[1:4])
  
  # remove clipping of x axis labels
  g2 <- ggplot_gtable(ggplot_build(g1))
  g2$layout$clip[g2$layout$name == "panel"] <- "off"
  return(arrangeGrob(g2))
}


# make single plot
makeepgplot=function(n,mn,prh1,prh2,prh3,ytitle){
  
  # format estimates
  # FECR was calculated as RR-1, but since effect is assumed
  # to be protective, switching to 1-RR
  # prh1=-prh1
  # prh2=-prh2
  # prh3=-prh3
  n=paste("(N=",format(n,big.mark=","),")",sep="")
  colnames(prh1)=c("rr","ub","lb","p-value")
  colnames(prh2)=c("rr","ub","lb","p-value")
  colnames(prh3)=c("rr","ub","lb","p-value")
  rr=apply(as.matrix(prh1),1,pt.est.ci.f,decimals=2,scale=1)
  rrh1=c("ref",rr)
  rrh2w=c("","ref","","",pt.est.ci.f(prh2[1,],decimals=2,scale=1),"","")
  rrh2s=c("","","ref","",pt.est.ci.f(prh2[2,],decimals=2,scale=1),"","")
  rrh2h=c("","","","ref",pt.est.ci.f(prh2[3,],decimals=2,scale=1),"","")
  rrh3wsh=c("","","","","ref","",pt.est.ci.f(prh3[1,],decimals=2,scale=1))
  rrh3n=c("","","","","","ref",pt.est.ci.f(prh3[2,],decimals=2,scale=1))
  cliney=as.numeric(mn$mn.f[1])
  
  # define color palette
  black = "#000004FF"
  blue = "#3366AA"
  teal = "#11AA99"
  green = "#66AA55"
  chartr = "#CCCC55"
  magent = "#992288"
  red = "#EE3333"
  orange = "#EEA722"
  yellow = "#FFEE33"
  grey = "#777777"
  cols=c(black,blue,teal,green,orange,red,magent,blue)
  
  colh2w=c("#5e5f60",cols[2],"#5e5f60","#5e5f60","#5e5f60","#5e5f60","#5e5f60")
  colh2s=c("#5e5f60","#5e5f60",cols[3],"#5e5f60","#5e5f60","#5e5f60","#5e5f60")
  colh2h=c("#5e5f60","#5e5f60","#5e5f60",cols[4],"#5e5f60","#5e5f60","#5e5f60")
  colh3wsh=c("#5e5f60","#5e5f60","#5e5f60","#5e5f60",cols[5],"#5e5f60","#5e5f60")
  colh3n=c("#5e5f60","#5e5f60","#5e5f60","#5e5f60","#5e5f60",cols[6],"#5e5f60")
  
  linediff=0.06
  secdiff=0.1
  
  if(ytitle=="Ascaris"){
    ytitlex=-0.73
    pbreaks=seq(0,10,1)
    plabels=seq(0,10,1)
    pylim=10
    pymin=0
    cmny2=cliney+((pylim-cliney)/20)*3
    cmny1=cliney+((pylim-cliney)/20)*1
  }
  if(ytitle=="Hookworm"){
    ytitlex=-0.58
    pbreaks=seq(0,1,.1)
    plabels=seq(0,1,.1)
    pylim=1
    pymin=0
    cmny2=cliney+((pylim-cliney)/20)*3.5
    cmny1=cliney+((pylim-cliney)/20)*1.5
  }
  if(ytitle=="Trichuris"){
    ytitlex=-0.7
    pbreaks=seq(0,1,.1)
    plabels=seq(0,1,.1)
    pylim=1
    pymin=0
    cmny2=cliney+((pylim-cliney)/20)*3.5
    cmny1=cliney+((pylim-cliney)/20)*1.5
  }
  
  g1=ggplot(mn,aes(x=x,y=mn.f))+
    geom_point(aes(col=x),alpha=0.7,size=1.5,show.legend=FALSE)+
    geom_errorbar(aes(ymin=lb,ymax=ub,col=x),width=0.11,show.legend=FALSE)+
    geom_hline(yintercept=cliney,linetype="dashed")+
    coord_cartesian(ylim = c(pymin, pylim),xlim=c(1,7)) +
    scale_color_manual("",values=cols)+
    ylab("Geometric mean\neggs per gram\nat 2-year\nfollow-up")+xlab("")+
    scale_y_continuous(breaks=pbreaks,labels=plabels)+
    theme_complete_bw() +
    theme(plot.margin = unit(c(7, 0.5, 0.5, 3.5), "lines"))+
    
    annotate(geom="text",x=0.65,y=cmny2,label="Control",size=2.5)+
    annotate(geom="text",x=0.62,y=cmny1,label="mean",size=2.5)+
    annotate(geom="text",x=ytitlex,y=pylim*(1.1+linediff*5+secdiff*2+0.15),label=ytitle,size=5.5)+
    annotate(geom="text",x=seq(1,7,1),y=pylim*(1.1+linediff*4+secdiff*2+0.16),
             label=c("Control","Water","Sanitation","Handwashing",
                     "Combined","Nutrition","Combined"),
             color=cols[1:7],size=3)+
    annotate(geom="text",x=seq(1,7,1),y=pylim*(1.1+linediff*4+secdiff*2+0.1),
             label=c("","","","","WSH","","Nutrition+WSH"),
             color=cols[1:7],size=3)+
    annotate(geom="text",x=seq(1,7,1),y=pylim*(1.1+linediff*4+secdiff*2+0.02),
             label=n,color="#5e5f60",size=3)+
    annotate(geom="text",x=-1.06,y=pylim*(1.1+linediff*5+secdiff*2+0.04), hjust=0, label="Fecal Egg Count Reduction",size=3.5)+
    annotate(geom="text",x=-1.06,y=pylim*(1.1+linediff*4+secdiff*2+0.02), hjust=0, label="(95% CI)",size=3.5)+
    annotate(geom="text",x=-0.1,y=pylim*(1.1+linediff*3+secdiff*2), label="Intervention vs. Control",size=2.7,col="#5e5f60")+
    
    annotate(geom="text",x=seq(1,7,1),y=pylim*(1.1+linediff*3+secdiff*2),label=rrh1,size=2.7,col="#5e5f60")+
    
    annotate(geom="text",x=0.22,y=pylim*(1.1+linediff*3+secdiff), label="WSH vs. W",size=2.7,col="#5e5f60")+
    annotate(geom="text",x=0.22,y=pylim*(1.1+linediff*2+secdiff), label="WSH vs. S",size=2.7,col="#5e5f60")+
    annotate(geom="text",x=0.22,y=pylim*(1.1+linediff+secdiff), label="WSH vs. H",size=2.7,col="#5e5f60")+
    annotate(geom="text",x=seq(1,7,1),y=pylim*(1.1+linediff*3+secdiff),label=rrh2w,size=2.7,col=colh2w)+
    annotate(geom="text",x=seq(1,7,1),y=pylim*(1.1+linediff*2+secdiff),label=rrh2s,size=2.7,col=colh2s)+
    annotate(geom="text",x=seq(1,7,1),y=pylim*(1.1+linediff+secdiff),label=rrh2h,size=2.7,col=colh2h)+
    
    annotate(geom="text",x=-.16,y=pylim*(1.1+linediff), label="Nutrition + WSH vs. WSH",size=2.7,col="#5e5f60")+
    annotate(geom="text",x=-.23,y=pylim*1.1, label="Nutrition + WSH vs. Nutrition",size=2.7,col="#5e5f60")+
    annotate(geom="text",x=seq(1,7,1),y=pylim*(1.1+linediff),label=rrh3wsh,size=2.7,col=colh3wsh)+
    annotate(geom="text",x=seq(1,7,1),y=pylim*1.1,label=rrh3n,size=2.7,col=colh3n)+
    annotate(geom="text",x=seq(1,7,1)+.15,y=mn$mn.f,label=mn$mn.f,size=2.7,col=cols[1:7])
  
  # remove clipping of x axis labels
  g2 <- ggplot_gtable(ggplot_build(g1))
  g2$layout$clip[g2$layout$name == "panel"] <- "off"
  
  return(arrangeGrob(g2))
}



# make single plot
makeemplot=function(prev,pr,ytitle,ylim,em1lab,em0lab,em1labs,em0labs){
  
  # format estimates
  rr1=apply(as.matrix(pr[pr$level==1,]),1,pt.est.ci.f,decimals=2,scale=1)
  rr1=c("ref",rr1)
  rr0=apply(as.matrix(pr[pr$level==0,]),1,pt.est.ci.f,decimals=2,scale=1)
  rr0=c("ref",rr0)
  
  prev$col=as.factor(seq(1,14,1))
  
  # define color palette
  black = "#000004FF"
  blue = "#3366AA"
  teal = "#11AA99"
  green = "#66AA55"
  orange = "#EEA722"
  red = "#EE3333"
  magent = "#992288"
  black2 = "#8B8B8C"
  blue2 = "#4284DB"
  teal2 = "#65DBCE"
  green2 = "#99DB88"
  orange2 = "#FAD796"
  red2 = "#F29494"
  magent2 = "#E381D5"
  cols=c(black,blue,teal,green,orange,red,magent)
  
  if(ytitle=="Ascaris"){
    ytitlex=-0.65
  }
  if(ytitle=="Hookworm"){
    ytitlex=-0.5
  }
  if(ytitle=="Trichuris"){
    ytitlex=-0.6
  }
  if(ytitle=="Any STH"){
    ytitlex=-0.59
  }
  
  pylim=ylim[2]
  linediff=0.06
  secdiff=0.1
  
  cliney=min(as.numeric(prev$prev.f[prev$tr=="Control"]))
  emlabmin=cliney*0.5
  
  em1=prev[prev$tr=="Control" & prev$level==1,"Prev"]
  em0=prev[prev$tr=="Control" & prev$level==0,"Prev"]
  
  g1=ggplot(prev,aes(x=x,y=prev.f,group=level))+
    geom_point(aes(col=x,shape=level),size=1.5,
               show.legend=FALSE,position=position_dodge(width=0.4))+
    geom_errorbar(aes(ymin=lb,ymax=ub,col=x),width=0.11,show.legend=FALSE,
                  position=position_dodge(width=0.4))+
    scale_shape_manual(name="Effect modifier level",values=c(21,19))+
    coord_cartesian(ylim = c(ylim[1], pylim),xlim=c(1,7))+
    scale_color_manual("",values=cols)+
    ylab("Prevalence\nat 2-year\nfollow-up (%)")+xlab("")+
    scale_y_continuous(breaks=seq(ylim[1],pylim,0.1),labels=seq(ylim[1],pylim,0.1))+
    theme_complete_bw() +
    theme(plot.margin = unit(c(5.5, 0.5, 0.5, 3.5), "lines"))+
    
    annotate(geom="text",x=-1.3,y=pylim*(1.1+linediff*3+secdiff*1+0.05),hjust=0,label=ytitle,size=5.5)+
    annotate(geom="text",x=seq(1,7,1),y=pylim*(1.1+linediff*2+secdiff),
             label=c("Control","Water","Sanitation","Handwashing",
                     "Combined","Nutrition","Combined"),
             color=cols[1:7],size=3)+
    annotate(geom="text",x=seq(1,7,1),y=pylim*(1.1+linediff*1+secdiff),
             label=c("","","","","WSH","","Nutrition+WSH"),
             color=cols[1:7],size=3)+
    annotate(geom="text",x=-1.3,y=pylim*(1.1+linediff*2+secdiff),hjust=0,  label="Prevalence Ratio (95% CI)",size=3.5)+
    annotate(geom="text",x=-1.3,y=pylim*(1.1+linediff*2),hjust=0, label="Intervention vs. Control",size=2.7)+
    annotate(geom="text",x=-1.3,y=pylim*(1.1+linediff*1), hjust=0,label=em0lab,size=2.7,col="#5e5f60")+
    annotate(geom="text",x=-1.3,y=pylim*(1.1), hjust=0,label=em1lab,size=2.7,col="#5e5f60")+
    
    annotate(geom="text",x=seq(1,7,1),y=pylim*(1.1+linediff),label=rr0,size=2.7,col="#5e5f60")+
    annotate(geom="text",x=seq(1,7,1),y=pylim*(1.1),label=rr1,size=2.7,col="#5e5f60")+
    
    annotate(geom="text",x=0.95,y=em0, hjust=0,label=em0labs,size=2.7,col="#5e5f60")+
    annotate(geom="text",x=1.15,y=em1, hjust=0,label=em1labs,size=2.7,col="#5e5f60")
  
  # remove clipping of x axis labels
  g2 <- ggplot_gtable(ggplot_build(g1))
  g2$layout$clip[g2$layout$name == "panel"] <- "off"
  # pdf("~/Box Sync/WASHB Parasites/Results/Figures/test.pdf",width=9,height=3.75)
  # grid.draw(g2)
  # dev.off()
  return(arrangeGrob(g2))
}




# make single plot
makeepgemplot=function(mn,pr,ytitle,ylim,em1lab,em0lab,em1labs,em0labs){
  
  mn$mn.f=as.numeric(mn$mn.f)
  
  # format estimates
  # FECR was calculated as RR-1, but since effect is assumed
  # to be protective, switching to 1-RR
  # pr$prf=pr$rr
  # pr$prf[pr$rr<0]=-pr$rr[pr$rr<0]
  # pr$lbf=pr$lb
  # pr$lbf[pr$rr<0]=-(pr$ub[pr$rr<0])
  # pr$ubf=pr$ub
  # pr$ubf[pr$rr<0]=-(pr$lb[pr$rr<0])
  
  # prf1=pr[pr$level==1,c("prf","lbf","ubf")]
  # colnames(prf1)=c("rr","lb","ub")
  # prf0=pr[pr$level==0,c("prf","lbf","ubf")]
  # colnames(prf0)=c("rr","lb","ub")
  
  prf1=pr[pr$level==1,c("rr","lb","ub")]
  prf0=pr[pr$level==0,c("rr","lb","ub")]
  
  # rr1=apply(as.matrix(prf1),1,pt.est.ci.f,decimals=2,scale=1)
  # rr1=c("ref",rr1)
  # rr0=apply(as.matrix(prf0),1,pt.est.ci.f,decimals=2,scale=1)
  # rr0=c("ref",rr0)
  
  rr1=apply(as.matrix(prf1),1,pt.est.ci.f,decimals=2,scale=1)
  rr1=c("ref",rr1)
  rr0=apply(as.matrix(prf0),1,pt.est.ci.f,decimals=2,scale=1)
  rr0=c("ref",rr0)
  
  # define color palette
  black = "#000004FF"
  blue = "#3366AA"
  teal = "#11AA99"
  green = "#66AA55"
  orange = "#EEA722"
  red = "#EE3333"
  magent = "#992288"
  black2 = "#8B8B8C"
  blue2 = "#4284DB"
  teal2 = "#65DBCE"
  green2 = "#99DB88"
  orange2 = "#FAD796"
  red2 = "#F29494"
  magent2 = "#E381D5"
  cols=c(black,blue,teal,green,orange,red,magent)
  
  linediff=0.06
  secdiff=0.1
  
  if(ytitle=="Ascaris"){
    pbreaks=seq(0,15,3)
    plabels=seq(0,15,3)
    pylim=15
    pymin=0
    ytitlex=-0.75
  }
  if(ytitle=="Hookworm"){
    ytitlex=-0.6
    pbreaks=seq(0,3,.5)
    plabels=seq(0,3,.5)
    pylim=3
    pymin=0
  }
  if(ytitle=="Trichuris"){
    ytitlex=-0.7
    pbreaks=seq(0,3,.5)
    plabels=seq(0,3,.5)
    pylim=3
    pymin=0
  }
  
  linediff=0.06
  secdiff=0.1
  
  cliney=min(as.numeric(mn$mn.f[mn$tr=="Control"]))
  emlabmin=cliney*0.5
  
  em1=mn[mn$tr=="Control" & mn$level==1,"geomean"]
  em0=mn[mn$tr=="Control" & mn$level==0,"geomean"]
  
  g1=ggplot(mn,aes(x=x,y=mn.f,group=level))+
    geom_point(aes(col=x,shape=level),size=1.5,
               show.legend=FALSE,position=position_dodge(width=0.4))+
    geom_errorbar(aes(ymin=lb,ymax=ub,col=x),width=0.11,show.legend=FALSE,
                  position=position_dodge(width=0.4))+
    scale_shape_manual(name="Effect modifier level",values=c(21,19))+
    coord_cartesian(ylim = c(pymin, pylim),xlim=c(1,7))+
    scale_color_manual("",values=cols)+
    ylab("Geometric mean\negg count\nat 2-year\nfollow-up (%)")+xlab("")+
    scale_y_continuous(breaks=pbreaks,labels=plabels)+
    theme_complete_bw() +
    theme(plot.margin = unit(c(5.5, 0.5, 0.5, 4), "lines"))+
    
    annotate(geom="text",x=-1.3,y=pylim*(1.1+linediff*3+secdiff*1+0.07),hjust=0,label=ytitle,size=5.5)+
    annotate(geom="text",x=seq(1,7,1),y=pylim*(1.1+linediff*2+secdiff),
             label=c("Control","Water","Sanitation","Handwashing",
                     "Combined","Nutrition","Combined"),
             color=cols[1:7],size=3)+
    annotate(geom="text",x=seq(1,7,1),y=pylim*(1.1+linediff*1+secdiff),
             label=c("","","","","WSH","","Nutrition+WSH"),
             color=cols[1:7],size=3)+
    annotate(geom="text",x=-1.3,y=pylim*(1.1+linediff*4.2),hjust=0, label="Adjusted Fecal Egg Count",size=3.5)+
    annotate(geom="text",x=-1.3,y=pylim*(1.1+linediff*3),hjust=0, label="Reduction [RR-1] (95% CI)",size=3.5)+   
    annotate(geom="text",x=-1.3,y=pylim*(1.1+linediff*2),hjust=0, label="Intervention vs. Control",size=2.7)+
    annotate(geom="text",x=-1.3,y=pylim*(1.1+linediff*1), hjust=0,label=em0lab,size=2.7,col="#5e5f60")+
    annotate(geom="text",x=-1.3,y=pylim*(1.1), hjust=0,label=em1lab,size=2.7,col="#5e5f60")+
    
    annotate(geom="text",x=seq(1,7,1),y=pylim*(1.1+linediff),label=rr0,size=2.7,col="#5e5f60")+
    annotate(geom="text",x=seq(1,7,1),y=pylim*(1.1),label=rr1,size=2.7,col="#5e5f60")+
    
    annotate(geom="text",x=0.95,y=em0, hjust=0,label=em0labs,size=2.7,col="#5e5f60")+
    annotate(geom="text",x=1.15,y=em1, hjust=0,label=em1labs,size=2.7,col="#5e5f60")
  
  
  # remove clipping of x axis labels
  g2 <- ggplot_gtable(ggplot_build(g1))
  g2$layout$clip[g2$layout$name == "panel"] <- "off"
  # pdf("~/Box Sync/WASHB Parasites/Results/Figures/test.pdf",width=9,height=3.75)
  # grid.draw(g2)
  # dev.off()
  return(arrangeGrob(g2))
}




# save 4 binary sth plots in one pdf
saveplot=function(alplot,hwplot,ttplot,sthplot,lab,res.dir){
  pdf(paste(res.dir,lab,".pdf",sep=""),width=9,height=15.5)
  grid.arrange(alplot,hwplot,ttplot,sthplot,ncol=1,nrow=4)
  dev.off()
}


# wrapper function
sth.bin.plot=function(aln,hwn,ttn,sthn,alprev,hwprev,ttprev,sthprev,alprh1,alprh2,hwprh1,hwprh2,ttprh1,ttprh2,sthprh1,sthprh2,lab,res.dir){
  # prep data for plotting
  alprev=sth.plot.prep(alprev)
  hwprev=sth.plot.prep(hwprev)
  ttprev=sth.plot.prep(ttprev)
  sthprev=sth.plot.prep(sthprev)
  
  alplot=makeplot(aln,alprev,alprh1,alprh2,ytitle="Ascaris",ylim)
  hwplot=makeplot(hwn,hwprev,hwprh1,hwprh2,ytitle="Hookworm",ylim)
  ttplot=makeplot(ttn,ttprev,ttprh1,ttprh2,ytitle="Trichuris",ylim)
  sthplot=makeplot(sthn,sthprev,sthprh1,sthprh2,ytitle="Any STH",ylim)
  
  saveplot(alplot,hwplot,ttplot,sthplot,lab)
}



# wrapper function
sth.epg.plot=function(aln,hwn,ttn,almn,hwmn,ttmn,algeoh1,algeoh2,algeoh3,hwgeoh1,hwgeoh2,hwgeoh3,ttgeoh1,ttgeoh2,ttgeoh3,lab){
  # prep data for plotting
  almn=sth.epg.plot.prep(almn)
  hwmn=sth.epg.plot.prep(hwmn)
  ttmn=sth.epg.plot.prep(ttmn)
  
  alplot=makeepgplot(aln,almn,algeoh1,algeoh2,algeoh3,ytitle="Ascaris")
  hwplot=makeepgplot(hwn,hwmn,hwgeoh1,hwgeoh2,hwgeoh3,ytitle="Hookworm")
  ttplot=makeepgplot(ttn,ttmn,ttgeoh1,ttgeoh2,ttgeoh3,ytitle="Trichuris")
  
  saveepgplot(alplot,hwplot,ttplot,lab)
}


# wrapper function
sth.epg.em.plot=function(almn1,hwmn1,ttmn1,almn0,hwmn0,ttmn0,alprh11,hwprh11,ttprh11,alprh10,hwprh10,ttprh10,figlab,em1lab,em0lab,em1labs,em0labs,ylim){
  almn1=sth.epg.em.plot.prep(almn1)
  hwmn1=sth.epg.em.plot.prep(hwmn1)
  ttmn1=sth.epg.em.plot.prep(ttmn1)
  
  almn0=sth.epg.em.plot.prep(almn0)
  hwmn0=sth.epg.em.plot.prep(hwmn0)
  ttmn0=sth.epg.em.plot.prep(ttmn0)
  
  almn=rbind(almn1,almn0)
  hwmn=rbind(hwmn1,hwmn0)
  ttmn=rbind(ttmn1,ttmn0)
  
  almn$level=as.factor(c(rep(1,7),rep(0,7)))
  hwmn$level=as.factor(c(rep(1,7),rep(0,7)))
  ttmn$level=as.factor(c(rep(1,7),rep(0,7)))
  
  alpr=rbind(alprh11,alprh10)
  hwpr=rbind(hwprh11,hwprh10)
  ttpr=rbind(ttprh11,ttprh10)
  
  alpr$level=c(rep(1,6),rep(0,6))
  hwpr$level=c(rep(1,6),rep(0,6))
  ttpr$level=c(rep(1,6),rep(0,6))
  
  alplot=makeepgemplot(almn,alpr,ytitle="Ascaris",em1lab=em1lab,em0lab=em0lab,em1labs=em1labs,em0labs=em0labs,ylim=ylim)
  hwplot=makeepgemplot(hwmn,hwpr,ytitle="Hookworm",em1lab=em1lab,em0lab=em0lab,em1labs=em1labs,em0labs=em0labs,ylim=ylim)
  ttplot=makeepgemplot(ttmn,ttpr,ytitle="Trichuris",em1lab=em1lab,em0lab=em0lab,em1labs=em1labs,em0labs=em0labs,ylim=ylim)
  
  figlab=paste("em-",figlab,sep="")
  saveepgplot(alplot,hwplot,ttplot,figlab)
}


# wrapper function
sth.bin.em.plot=function(alprev1,hwprev1,ttprev1,sthprev1,alprev0,hwprev0,ttprev0,sthprev0,alprh11,hwprh11,ttprh11,sthprh11,alprh10,hwprh10,ttprh10,sthprh10,figlab,em1lab,em0lab,em1labs,em0labs,ylim){
  alprev1=sth.plot.prep(alprev1)
  hwprev1=sth.plot.prep(hwprev1)
  ttprev1=sth.plot.prep(ttprev1)
  sthprev1=sth.plot.prep(sthprev1)
  
  alprev0=sth.plot.prep(alprev0)
  hwprev0=sth.plot.prep(hwprev0)
  ttprev0=sth.plot.prep(ttprev0)
  sthprev0=sth.plot.prep(sthprev0)
  
  alprev=rbind(alprev1,alprev0)
  hwprev=rbind(hwprev1,hwprev0)
  ttprev=rbind(ttprev1,ttprev0)
  sthprev=rbind(sthprev1,sthprev0)
  
  alprev$level=as.factor(c(rep(1,7),rep(0,7)))
  hwprev$level=as.factor(c(rep(1,7),rep(0,7)))
  ttprev$level=as.factor(c(rep(1,7),rep(0,7)))
  sthprev$level=as.factor(c(rep(1,7),rep(0,7)))
  
  alpr=rbind(alprh11,alprh10)
  hwpr=rbind(hwprh11,hwprh10)
  ttpr=rbind(ttprh11,ttprh10)
  sthpr=rbind(sthprh11,sthprh10)
  
  alpr$level=c(rep(1,6),rep(0,6))
  hwpr$level=c(rep(1,6),rep(0,6))
  ttpr$level=c(rep(1,6),rep(0,6))
  sthpr$level=c(rep(1,6),rep(0,6))
  
  alplot=makeemplot(alprev,alpr,ytitle="Ascaris",em1lab=em1lab,em0lab=em0lab,em1labs=em1labs,em0labs=em0labs,ylim=ylim)
  hwplot=makeemplot(hwprev,hwpr,ytitle="Hookworm",em1lab=em1lab,em0lab=em0lab,em1labs=em1labs,em0labs=em0labs,ylim=ylim)
  ttplot=makeemplot(ttprev,ttpr,ytitle="Trichuris",em1lab=em1lab,em0lab=em0lab,em1labs=em1labs,em0labs=em0labs,ylim=ylim)
  sthplot=makeemplot(sthprev,sthpr,ytitle="Any STH",em1lab=em1lab,em0lab=em0lab,em1labs=em1labs,em0labs=em0labs,ylim=ylim)
  
  figlab=paste("em-",figlab,sep="")
  saveplot(alplot,hwplot,ttplot,sthplot,figlab)
}