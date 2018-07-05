#######################################
# WASH Benefits Bangladesh STH KK qPCR validation
# inputs for bayesian lca
# simulation with known parameters
#######################################
rm(list=ls())
library(dplyr)
library(nimble)
library(reshape2)
library(ggplot2)

# set.seed(123)
# se1=rbeta(1000,shape1=13.3, shape2=2)
# hist(se1)
# se2=rbeta(1000,shape1=40, shape2=40)
# hist(se2)
# 
# Mode <- function(x) {
#   ux <- unique(x)
#   ux[which.max(tabulate(match(x, ux)))]
# }
# 
# Mode(se)
# quantile(se,c(0,0.05,1))
# 
# 
set.seed(123)
sp=rbeta(1000,shape1=21.2, shape2=1)
hist(sp)

set.seed(123)
sp=rbeta(1000,shape1=5, shape2=5)
hist(sp)


set.seed(123)
sp=rbeta(1000,shape1=5, shape2=20)
hist(sp)
# 
# 
# Mode(sp)
# quantile(sp,c(0,0.05,1))
sp=rbeta(1000,shape1=21.2, shape2=1)
hist(sp)

# Easton: 95 % confidence that the speci- ficity of each method was greater than 95 %
# uniform distribution (does better)
set.seed(123)
d=data.frame(sp=runif(1000,min=0.9475, max=1))
d$over95=ifelse(d$sp>0.95,1,0)
prop.table(table(d$over95))
hist(d$sp)

# gaussian distribution
set.seed(123)
d=data.frame(sp=rnorm(1000,mean=0.975, sd=0.0152))
d$over95=ifelse(d$sp>0.95,1,0)
prop.table(table(d$over95))
hist(d$sp)

# Easton: 95 % con- fidence that the sensitivity of qPCR was greater than 80 %.”
# uniform distribution (does better)
set.seed(123)
d=data.frame(se=runif(1000,min=0.7895, max=1))
d$over95=ifelse(d$se>0.8,1,0)
prop.table(table(d$over95))
hist(d$se)

# gaussian distribution
set.seed(123)
d=data.frame(se=rnorm(1000,mean=0.9, sd=0.1255))
d$over95=ifelse(d$se>0.8,1,0)
prop.table(table(d$over95))
hist(d$se)





# assume sensitivity = 0.95, sensitivity = 0.98
# true prevalence = 0.3, correlation
# use similar n

sens.q=0.95
spec.q=0.98
sens.k=0.5
spec.k=0.98
pi=0.25
n=2800

qplus=n*((pi*sens.q)+(1-pi)*(1-spec.q))
qminus=n*((pi*(1-sens.q))+((1-pi)*spec.q))

d <- data.frame(matrix(ncol = 1, nrow = n))
colnames(d)="q"

set.seed(123)
d <- d %>%
  mutate(y=rbinom(n,1,pi),
         q=rbinom(n,1,y*sens.q+(1-y)*(1-spec.q)),
         k=rbinom(n,1,y*sens.k+(1-y)*(1-spec.k)))

# observed values
#prevalence
mean(d$y)

# sensitivity qpcr
prop.table(table(d$q[d$y==1]))[2]
# specificity qpcr
prop.table(table(d$q[d$y==0]))[1]

# sensitivity kk
prop.table(table(d$k[d$y==1]))[2]
# specificity kk
prop.table(table(d$k[d$y==0]))[1]

#covariance
cov(d$k[d$y==1],d$q[d$y==1])
cov(d$k[d$y==0],d$q[d$y==0])

# correlation 
cor(d$k[d$y==1],d$q[d$y==1])
cor(d$k[d$y==0],d$q[d$y==0])

x11=table(d$q,d$k)[2,2]
x10=table(d$q,d$k)[1,2]
x01=table(d$q,d$k)[2,1]
x00=table(d$q,d$k)[1,1]


sim.data=list(n=n, x=c(x11, x10, x01, x00))
sim.data
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
  Se1 ~ dbeta(1,1) 
  Sp1 ~ dbeta(1,1) 
  Se2 ~ dbeta(1,1) 
  Sp2 ~ dbeta(1,1) 
  
  covDn ~ dunif(0, uc)
  covDp ~ dunif(0, us)
  rhoD <- covDp / sqrt(Se1*(1-Se1)*Se2*(1-Se2))
  rhoDc <- covDn / sqrt(Sp1*(1-Sp1)*Sp2*(1-Sp2))
})

# define initial values
set.seed(123)
hw.inits=list(pi=0.22, Se1=0.526, Se2=0.98, Sp1=0.986, Sp2=0.97, 
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
lca.model.hw <- nimbleModel(lca.hw, inits=hw.inits, data=sim.data)

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
thin=100
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
hw.inits=list(list(pi=0.22, Se1=0.526, Se2=0.98, Sp1=0.986, Sp2=0.97,
                   rhoD=0.5, rhoDc=0.5,
                   p=c(0.526*0.22+0.98*0.22,
                       0.526*0.22+0.97*(1-0.22),
                       0.986*(1-0.22)+0.98*0.22,
                       0.986*(1-0.22)+0.97*(1-0.22))),
              list(pi=0.5, Se1=0.5, Se2=0.5, Sp1=0.5, Sp2=0.5,
                   rhoD=0.5, rhoDc=0.5,
                   p=c(0.5,0.5,0.5,0.5)),
              list(pi=0.3, Se1=0.3, Se2=0.3, Sp1=0.3, Sp2=0.3,
                   rhoD=0.3, rhoDc=0.3,
                   p=c(0.3,0.3,0.3,0.3)),
              list(pi=0.4, Se1=0.4, Se2=0.4, Sp1=0.4, Sp2=0.4,
                   rhoD=0.4, rhoDc=0.4,
                   p=c(0.4,0.4,0.4,0.4)))
set.seed(12345)
mcmc.hw.out <- runMCMC(C.hw.MCMC, niter = niter, nchains = 4, inits=hw.inits)

names(mcmc.hw.out)

#--------------------------------------
# Plot features of posterior distribution
#--------------------------------------
samples.hw.f=as.data.frame(do.call(rbind,mcmc.hw.out))
samples.hw.f$chain=as.factor(c(rep(1,niter/thin),
                               rep(2,niter/thin),rep(3,niter/thin),
                               rep(4,niter/thin)))

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
