library(phylofactor)
source('')
setwd('Mammals/')
pf <- readRDS('pf_BodyMass')

# 
# prob.fcn <- function(obj,lambda) obj^lambda/(sum(obj^lambda))
# obj <- sapply(getPhyloGroups(tree),TestFunction,tree,Z) %>% sort
# prob.fcn(obj,-6)[1]  ## 1/4 chance of drawing the winning edge.
# 

binframe <- function(pff,factor=pff$nfactors) data.frame('species'=pff$tree$tip.label,'Z'=pff$Data,'Bin'=factor(phylobin(bins(pff$basis[,1:factor,drop=F]))))
getDev <- function(pff,factor=pff$nfactors){
  return(anova(glm(Z~Bin,data=binframe(pff,factor)))['Bin','Deviance'])
}

reps=3000
pf.max=NULL
amax=0
abase <- getDev(pf,5)
abasew <- getDev(pf.W,5)
for (rr in 1:reps){
  pff <- twoSampleFactor(Z,tree,nfactors=5,Metropolis = T,lambda=6,ncores=7)
  a <- getDev(pff)
  if (a>amax){
    pf.max <- pff
    amax <- a
  }
}

getDev(pf)
getDev(pf.W)
getDev(pf.max)

B <- binframe(pf,5)
Bmax <- binframe(pf,5)
anova(glm(Z~Bin,data=B))

nn=5
s <- summary(pf.max,Tax,nn)
s$taxa.split$group1
mean(Z[pf.max$groups[[nn]][[1]]]) %>% exp
sapply(pf$bins,FUN=function(ix,Z) exp(mean(Z[ix])),Z) %>% signif(2) %>% sort
# sapply(pff$bins,FUN=function(ix,Z) exp(mean(Z[ix])),Z) %>% signif(2) %>% sort

save(list=ls(),file='Mammal_body_mass_Metropolis_Hastings_workspace')

pf$tree$edge.length <- rep(1,Nedge(pf$tree))
pf.tree(pf,factor.map = data.frame('factors'=1:5,'group'=1))
