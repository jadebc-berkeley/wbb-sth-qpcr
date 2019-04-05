#######################################
# WASH Benefits Bangladesh STH KK qPCR validation
# venn diagram of co-infection
#######################################
rm(list=ls())

# configure directories, load libraries and base functions
source(paste0(here::here(), "/0-config.R"))

# load data 
load(paste0(data_dir,"qdata.RData"))

# get n's - KK
get_kk_n = function(data, al, hw, tt){
  x = data %>% filter(
    alkk == al,
    hwkk == hw, 
    ttkk == tt
  ) %>%
    summarise(n = n())
  return(as.numeric(x))
}

# get n's - qPCR
get_q_n = function(data, al, hw, tt){
  x = data %>% filter(
    positive.Al2 == al,
    positive.Hw == hw, 
    positive.Tt == tt
  ) %>%
    summarise(n = n())
  return(as.numeric(x))
}

# get n's - qPCR Hw
get_hw_n = function(data, ac, na, ad){
  x = data %>% filter(
    positive.Ac == ac,
    positive.Na == na, 
    positive.Ad == ad
  ) %>%
    summarise(n = n())
  return(as.numeric(x))
}


# ---------------------------------------
# kato katz plot
# ---------------------------------------
kk=qdata[,c("alkk","hwkk","ttkk")]
kk$coinf=rowSums(kk,na.rm=TRUE)

# intersections
alhw=get_kk_n(kk, al = 1, hw = 1, tt = 0)
hwtt=get_kk_n(kk, al = 0, hw = 1, tt = 1)
altt=get_kk_n(kk, al = 1, hw = 0, tt = 1)
alhwtt=get_kk_n(kk, al = 1, hw = 1, tt = 1)

altot=get_kk_n(kk, al = 1, hw = 0, tt = 0)
hwtot=get_kk_n(kk, al = 0, hw = 1, tt = 0)
tttot=get_kk_n(kk, al = 0, hw = 0, tt = 1)

# circles proportional to size 
v <- venneuler(c(Ascaris=altot, Hookworm=hwtot, Trichuris=tttot,
                 "Ascaris&Hookworm"=alhw,"Hookworm&Trichuris"=hwtt,
                 "Trichuris&Ascaris"=altt,"Ascaris&Hookworm&Trichuris"=alhwtt))
v$labels <- c(
  "","",""
)

cb.lightorange="#E69F00"
cb.blue= "#56B4E9"
cb.green="#009E73"
cb.orange="#D55E00"
cb.pink="#CC79A7"
cb.dblue="#005787"
purple="#9E4AED"

fontinflate = 1.2

pdf(file=paste0(fig_dir, "wbb-venn-kk.pdf"),
    width=6,height=6)
  plot(v, col=c(cb.lightorange, cb.blue, cb.green))
  # titles
  text(x=v$centers[1,1],y=v$centers[1,2]+0.32, labels=expression(paste(italic("A. lumbricoides"))),
       cex = fontinflate)
  text(x=v$centers[1,1]+0.3,y=v$centers[1,2]-0.28, labels="Hookworm",
       cex = fontinflate)
  text(x=v$centers[1,1]+0.31,y=v$centers[1,2]+0.22, labels=expression(paste(italic("T. trichiura"))),
       cex = fontinflate)
  
  # non intersections
  text(x=v$centers[1,1],y=v$centers[1,2], labels=paste(altot), cex = fontinflate)
  text(x=v$centers[2,1]+0.06,y=v$centers[2,2]-0.05, labels=paste(hwtot), cex = fontinflate)
  text(x=v$centers[3,1]+0.06,y=v$centers[3,2]+0.02, labels=paste(tttot), cex = fontinflate)
  
  # intersections
  text(x=v$centers[1,1]+0.23,y=v$centers[1,2]+0.07, labels=paste(altt), cex = fontinflate)
  text(x=v$centers[1,1]+0.18,y=v$centers[1,2]-0.13, labels=paste(alhw), cex = fontinflate)
  text(x=v$centers[1,1]+0.23,y=v$centers[1,2]-0.035, labels=paste(alhwtt), cex = fontinflate)
  text(x=v$centers[1,1]+0.32,y=v$centers[1,2]-0.047, labels=paste(hwtt), cex = fontinflate)
dev.off()

# ---------------------------------------
# qpcr plot
# ---------------------------------------
# counts
altot=get_q_n(data = qdata, al = 1, hw = 0, tt = 0)
hwtot=get_q_n(data = qdata, al = 0, hw = 1, tt = 0)
tttot=get_q_n(data = qdata, al = 0, hw = 0, tt = 1)

# intersections
alhw=get_q_n(data = qdata, al = 1, hw = 1, tt = 0)
altt=get_q_n(data = qdata, al = 1, hw = 0, tt = 1)
hwtt=get_q_n(data = qdata, al = 0, hw = 1, tt = 1)

alhwtt=get_q_n(data = qdata, al = 1, hw = 1, tt = 1)

qpcr <- venneuler(c(Ascaris=altot, Hookworm=hwtot, Trichuris=tttot,
                 "Ascaris&Hookworm"=alhw,"Hookworm&Trichuris"=hwtt,
                 "Trichuris&Ascaris"=altt,"Ascaris&Hookworm&Trichuris"=alhwtt))
qpcr$labels <- c("","","")

centers=qpcr$centers

pdf(file=paste0(fig_dir, "wbb-venn-qpcr.pdf"),
    width=6,height=6)
  plot(qpcr, col=c(cb.lightorange, cb.blue,cb.green))
  # titles
  text(x=centers[1,1]+0.0,y=centers[1,2]+0.26, labels=expression(paste(italic("A. lumbricoides"))), cex = 1.2)
  text(x=centers[2,1]-0.2,y=centers[2,2]-0.22, labels="Hookworm", cex = fontinflate)
  text(x=centers[3,1]-0.12,y=centers[3,2]+0.2, labels=expression(paste(italic("T. trichiura"))), cex = fontinflate)
  
  # non intersections
  text(x=centers[1,1]+0.04,y=centers[1,2]+0.08, labels=paste(altot), cex = fontinflate)
  text(x=centers[2,1]-0.01,y=centers[2,2]-0.1, labels=paste(hwtot), cex = fontinflate)
  text(x=centers[3,1]-0.08,y=centers[3,2]+0.02, labels=paste(tttot), cex = fontinflate)
  
  # intersections
  text(x=centers[1,1]-0.16,y=centers[1,2]+0.02, labels=paste(altt), cex = fontinflate)
  text(x=centers[1,1]-0.13,y=centers[1,2]-0.1, labels=paste(alhwtt), cex = fontinflate)
  text(x=centers[3,1]-0.02,y=centers[3,2]-0.1, labels=paste(hwtt), cex = fontinflate)
  text(x=centers[1,1]+0.02,y=centers[1,2]-0.15, labels=paste(alhw), cex = fontinflate) 

dev.off()


# ---------------------------------------
# qpcr hw plot
# ---------------------------------------
# counts
actot=get_hw_n(qdata, ac = 1, na = 0, ad = 0)
natot=get_hw_n(qdata, ac = 0, na = 1, ad = 0)
adtot=get_hw_n(qdata, ac = 0, na = 0, ad = 1)

# intersections
acna=get_hw_n(qdata, ac = 1, na = 1, ad = 0)
acad=get_hw_n(qdata, ac = 1, na = 0, ad = 1)
naad=get_hw_n(qdata, ac = 0, na = 1, ad = 1)

acnaad=get_hw_n(qdata, ac = 1, na = 1, ad = 1)

qpcr.hw <- venneuler(c("A. ceylanicum"=actot, 
                    "N. americanus"=natot, 
                    "A. duodenale"=adtot,
                    "A. ceylanicum&N. americanus"=acna,
                    "A. ceylanicum&A. duodenale"=acad,
                    "N. americanus&A. duodenale"=naad,
                    "A. ceylanicum&N. americanus&A. duodenale"=acnaad))

qpcr.hw$labels <- c("","","")

centers=qpcr.hw$centers

pdf(file=paste0(fig_dir,"wbb-venn-qpcr-hw.pdf"),
    width=6,height=6)
plot(qpcr.hw, col=c(cb.lightorange, cb.blue, purple,cb.green))
# titles
text(x=centers[1,1],y=centers[1,2]+0.15, labels=expression(paste(italic("A. ceylanicum"))), cex = fontinflate)
text(x=centers[2,1],y=centers[2,2]-0.315, labels=expression(paste(italic("N. americanus"))), cex = fontinflate)
text(x=centers[3,1]-0.08,y=centers[3,2]+0.04, labels=expression(paste(italic("A. duodenale"))), cex = fontinflate)
# non intersections
text(x=centers[1,1],y=centers[1,2]+0.03, labels=paste(actot), cex = fontinflate)
text(x=centers[2,1],y=centers[2,2]-0.02, labels=paste(natot), cex = fontinflate)
# intersections
text(x=centers[1,1]+0.0,y=centers[1,2]-0.09, labels=paste(acna), cex = fontinflate)
text(x=centers[1,1]-0.245,y=centers[1,2]-0.165, labels=paste(naad),cex=0.75)
dev.off()
