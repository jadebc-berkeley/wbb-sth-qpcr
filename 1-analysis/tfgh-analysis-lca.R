#######################################
# WASH Benefits Bangladesh STH KK qPCR validation
# bayesian lca analysis to estimate
# true prevalence, sensitivity, specificity
#######################################
rm(list=ls())
library(nimble)
library(reshape2)
library(ggplot2)
load("~/Dropbox/WASH-B-STH-Add-on/TFGH/Data/RData/qdata.RData")


# define data
adata=qdata[!is.na(qdata$positive.Al) & !is.na(qdata$alkk),]
al.x11=nrow(adata[adata$positive.Al==1 & adata$alkk==1,])
al.x10=nrow(adata[adata$positive.Al==0 & adata$alkk==1,])
al.x01=nrow(adata[adata$positive.Al==1 & adata$alkk==0,])
al.x00=nrow(adata[adata$positive.Al==0 & adata$alkk==0,])
al.n=sum(al.x11,al.x10,al.x01,al.x00)

al.data=list(n=al.n, x=c(al.x11, al.x10, al.x01, al.x00))

hdata=qdata[!is.na(qdata$positive.Hw) & !is.na(qdata$hwkk),]
hw.x11=nrow(hdata[hdata$positive.Hw==1 & hdata$hwkk==1,])
hw.x10=nrow(hdata[hdata$positive.Hw==0 & hdata$hwkk==1,])
hw.x01=nrow(hdata[hdata$positive.Hw==1 & hdata$hwkk==0,])
hw.x00=nrow(hdata[hdata$positive.Hw==0 & hdata$hwkk==0,])
hw.n=sum(hw.x11,hw.x10,hw.x01,hw.x00)

hw.data=list(n=hw.n, x=c(hw.x11, hw.x10, hw.x01, hw.x00))

tdata=qdata[!is.na(qdata$positive.Tt) & !is.na(qdata$ttkk),]
tt.x11=nrow(tdata[tdata$positive.Tt==1 & tdata$ttkk==1,])
tt.x10=nrow(tdata[tdata$positive.Tt==0 & tdata$ttkk==1,])
tt.x01=nrow(tdata[tdata$positive.Tt==1 & tdata$ttkk==0,])
tt.x00=nrow(tdata[tdata$positive.Tt==0 & tdata$ttkk==0,])
tt.n=sum(tt.x11,tt.x10,tt.x01,tt.x00)

tt.data=list(n=tt.n, x=c(tt.x11, tt.x10, tt.x01, tt.x00))

#--------------------------------------
# Build BUGS model
#--------------------------------------
lca.hw <- nimbleCode({
  # stochastic node
  x[1:4] ~ dmulti(p[1:4], n) 
  
  # deterministic nodes
  p[1] <- pi*(Se1*Se2+covDp) + (1-pi)*((1-Sp1)*(1-Sp2)+covDn)
  p[2] <- pi*(Se1*(1-Se2)-covDp) + (1-pi)*((1-Sp1)*Sp2-covDn)
  p[3] <- pi*((1-Se1)*Se2-covDp) + (1-pi)*(Sp1*(1-Sp2)-covDn)
  p[4] <- pi*((1-Se1)*(1-Se2)+covDp) + (1-pi)*(Sp1*Sp2+covDn)
  
  us <- min(Se1,Se2) - Se1*Se2
  uc <- min(Sp1,Sp2) - Sp1*Sp2
  
  # prior values
  pi ~ dbeta(1,1) 
  Se1 ~ dbeta(3,3) 
  Sp1 ~ dunif(0.9475,1) 
  Se2 ~ dunif(0.7895,1) 
  Sp2 ~ dunif(0.9475,1)  
  
  covDn ~ dunif(0, uc)
  covDp ~ dunif(0, us)
  rhoD <- covDp / sqrt(Se1*(1-Se1)*Se2*(1-Se2))
  rhoDc <- covDn / sqrt(Sp1*(1-Sp1)*Sp2*(1-Sp2))
})

# define initial values
set.seed(123)
hw.inits=list(pi=0.22, Se1=0.33, Se2=0.99, Sp1=0.986, Sp2=0.97, 
              rhoD=0.5, rhoDc=0.5, 
              p=c(0.526*0.22+0.98*0.22, 
                  0.526*0.22+0.97*(1-0.22), 
                  0.986*(1-0.22)+0.98*0.22, 
                  0.986*(1-0.22)+0.97*(1-0.22)), 
              covDp=runif(1,min=0,max=0.526-(0.526*0.98)),
              covDn=runif(1,min=0,max=0.97-(0.986*0.97)))


#--------------------------------------
# Configure model (nimbleModel)
#--------------------------------------
# process BUGS model code
lca.model.hw <- nimbleModel(lca.hw, inits=hw.inits, data=hw.data)

# list nodes without initial values
lca.model.hw$initializeInfo()

# display initial values
lca.model.hw[["p"]] 

# simulate values 
lca.model.hw$simulate("pi")

# list node names
lca.model.hw$getNodeNames()

# plot DAG
lca.model.hw$plotGraph()

# show all dependencies of pi terminating in stochastic nodes
lca.model.hw$getDependencies(c("pi"))

c.lca.model.hw=compileNimble(lca.model.hw)

#--------------------------------------
# Configure MCMC 
#--------------------------------------
monitors=c("pi","Se1","Se2","Sp1","Sp2","rhoD", "rhoDc")
thin=10
lca.hw.Conf <- configureMCMC(lca.model.hw,
      monitors=monitors,
      thin=thin, print = TRUE)
lca.hw.Conf$printSamplers()
# default is Metropolis-Hastings sampler
# how do I specify gibbs sampling? help(sampler) - see example of building with dmulti

# thin indicates when samples should be recorded- what should this be?

#--------------------------------------
# Build MCMC 
#--------------------------------------
hw.MCMC <- buildMCMC(lca.hw.Conf, enableWAIC = TRUE)

#--------------------------------------
# Compile MCMC 
#--------------------------------------
C.hw.MCMC <- compileNimble(hw.MCMC, project=lca.model.hw)

#--------------------------------------
# Run MCMC
#--------------------------------------
niter=100000

# specify initial values for each chain
hw.inits=list(list(pi=0.22, Se1=0.33, Se2=0.99, Sp1=0.986, Sp2=0.97, 
                   rhoD=0.5, rhoDc=0.5, 
                   p=c(0.526*0.22+0.98*0.22, 
                       0.526*0.22+0.97*(1-0.22), 
                       0.986*(1-0.22)+0.98*0.22, 
                       0.986*(1-0.22)+0.97*(1-0.22)), 
                   covDp=runif(1,min=0,max=0.526-(0.526*0.98)),
                   covDn=runif(1,min=0,max=0.97-(0.986*0.97))),
              
              list(pi=0.22, Se1=0.33, Se2=0.99, Sp1=0.986, Sp2=0.97, 
                   rhoD=0.5, rhoDc=0.5, 
                   p=c(0.526*0.22+0.98*0.22, 
                       0.526*0.22+0.97*(1-0.22), 
                       0.986*(1-0.22)+0.98*0.22, 
                       0.986*(1-0.22)+0.97*(1-0.22)), 
                   covDp=runif(1,min=0,max=0.526-(0.526*0.98)),
                   covDn=runif(1,min=0,max=0.97-(0.986*0.97))),
              
              list(pi=0.22, Se1=0.33, Se2=0.99, Sp1=0.986, Sp2=0.97, 
                   rhoD=0.5, rhoDc=0.5, 
                   p=c(0.526*0.22+0.98*0.22, 
                       0.526*0.22+0.97*(1-0.22), 
                       0.986*(1-0.22)+0.98*0.22, 
                       0.986*(1-0.22)+0.97*(1-0.22)), 
                   covDp=runif(1,min=0,max=0.526-(0.526*0.98)),
                   covDn=runif(1,min=0,max=0.97-(0.986*0.97))),
              
              list(pi=0.22, Se1=0.33, Se2=0.99, Sp1=0.986, Sp2=0.97, 
                   rhoD=0.5, rhoDc=0.5, 
                   p=c(0.526*0.22+0.98*0.22, 
                       0.526*0.22+0.97*(1-0.22), 
                       0.986*(1-0.22)+0.98*0.22, 
                       0.986*(1-0.22)+0.97*(1-0.22)), 
                   covDp=runif(1,min=0,max=0.526-(0.526*0.98)),
                   covDn=runif(1,min=0,max=0.97-(0.986*0.97))))
set.seed(12345)
mcmc.hw.out <- runMCMC(C.hw.MCMC, niter = niter, nchains = 4, inits=hw.inits)

names(mcmc.hw.out)

#--------------------------------------
# Get mean and percentiles of posterior distribution, combining chains
#--------------------------------------
samples.hw.f=as.data.frame(do.call(rbind,mcmc.hw.out))
samples.hw.f$chain=as.factor(c(rep(1,niter/thin),
                               rep(2,niter/thin),rep(3,niter/thin),
                               rep(4,niter/thin)))

hw.pi=c(mean(samples.hw.f$pi),quantile(samples.hw.f$pi,probs=c(0.025, 0.975)))
hw.Se1=c(mean(samples.hw.f$Se1),quantile(samples.hw.f$Se1,probs=c(0.025, 0.975)))
hw.Se2=c(mean(samples.hw.f$Se2),quantile(samples.hw.f$Se2,probs=c(0.025, 0.975)))
hw.Sp1=c(mean(samples.hw.f$Sp1),quantile(samples.hw.f$Sp1,probs=c(0.025, 0.975)))
hw.Sp2=c(mean(samples.hw.f$Sp2),quantile(samples.hw.f$Sp2,probs=c(0.025, 0.975)))
hw.rhoD=c(mean(samples.hw.f$rhoD),quantile(samples.hw.f$rhoD,probs=c(0.025, 0.975)))
hw.rhoDc=c(mean(samples.hw.f$rhoDc),quantile(samples.hw.f$rhoDc,probs=c(0.025, 0.975)))

hw.out=t(cbind(hw.pi,hw.Se1,hw.Sp1,hw.Sp1,hw.Sp2,hw.rhoD,hw.rhoDc))
hw.out=data.frame(monitors,hw.out)
colnames(hw.out)=c("monitors","mean","lb","ub")

#--------------------------------------
# Plot features of posterior distribution
#--------------------------------------
# convert samples to long format for plotting
samples.hw.l=melt(samples.hw.f,id.vars=c("chain"))

# plot history of each parameter
numMonitors=7
xseq=rep(rep(seq(1,niter/thin),numMonitors),4)

ggplot(samples.hw.l, aes(x=xseq,y=value,color=chain))+
  geom_line()+facet_wrap(~variable, scales="free",ncol=2)+theme_bw()+
  xlab("Iteration")+ylab("Parameter estimate")

# plot density of each parameter
ggplot(samples.hw.l, aes(x=value))+geom_density(aes(col=chain))+
  facet_wrap(~variable, scales="free")+theme_bw()+
  xlab("Kernal density")+xlab("Parameter estimate")

# plot autocorrelation for each parameter
par(mfrow = c(1, 5), mai = c(.6, .4, .1, .2))
acf(samples.hw.f[, "pi"])
acf(samples.hw.f[, "Se1"]) 
acf(samples.hw.f[, "Se2"]) 
acf(samples.hw.f[, "Sp1"]) 
acf(samples.hw.f[, "Sp2"]) 

# gelman-rubin diagnostic
# DIC
# MC error

