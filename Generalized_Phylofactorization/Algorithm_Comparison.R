#### In this script, we illustrate the pro's/con's of different algorithms
#### for glm phylofactorization 
library(phylofactor)
library(parallel)
library(tidyverse)
require(phangorn)
setwd('Generalized_Phylofactorization')
## Given coefficient matrix B, and likelihood L, the metrics of success will be:
## (1) ||B-B_est||  approximation of coefficient matrix
## (2) %B           % bins correct
## (3) L(k)         likelihood of model for factor k
## (4) L(tot,k)     global likelihood of total model for factors 1,...,k
## (5) time         simulation time.

## Algorithms compared will be:
## |v'B|      -      projection of B onto contrast vector
## |v's(B)|   -      projection of standardized B onto contrast vector
## L(phylo)   -      deviance of phylo*(PartitioningVariables)

##### Our formula will be
formula = Data~phylo*(y+z)
##### where species have species-specific mean abundance, const, and phylogenetic factors differentiate respons to y and z. 
##### and our partitioning variables will be y and z. Both will have delta = 0.3 for first factor, delta = 0.2 for the second.

##### We will have some multicolinearity between z and y in an effort to illustrate the tradeoffs between these algorithms. The multicolinarity is simulated via
##### z=y+rnorm(n)

##### species will have different mean abundances, "const", and there will be two edges with a "phylo" effect.

set.seed(1)
m=50 ## number of species
n=40 ## number of samples
tree <- rtree(m)

clade1 <- phangorn::Descendants(tree,75,'tips')[[1]]
clade2 <- phangorn::Descendants(tree,53,'tips')[[1]]
Grps <- getPhyloGroups(tree)
V <- sapply(Grps,ilrvec,m)
ilogit <- function(eta) 1/(1+exp(-eta))

const=rnorm(m)   ### logarithm of mean abundance, fixed across simulations

datasim <- function(n=40,m=50,delta1=0.3,delta2=0.3,cnst=const,c1=clade1,c2=clade2,tr=tree){
  y <- rnorm(n)
  z <- y+rnorm(n)
  etas <- t(matrix(rep(cnst,times=n),nrow=m,ncol=n,byrow=F))
  etas <- etas-.1*y-.1*z
  etas[,c1] <- etas[,c1]-delta1*y+2*delta1*z
  etas[,c2] <- etas[,c2]-delta2*y+2*delta2*z
  done=F
  while(!done){
    Data <- matrix(rbinom(m*n,size=1,prob=ilogit(etas)),nrow=m,ncol=n,byrow=T)
    rs <- rowSums(Data)
    if (min(rs)>0 & max(rs)<n){
      done=T
    }
  }
  colnames(Data) <- as.character(1:n)
  rownames(Data) <- tr$tip.label
  X <- data.frame(y,z)
  return(list('Data'=Data,'X'=X))
}

Bproj <- function(V,coefs,WinningPoolSize=1){
  omegas <- apply((coefs[c('y','z'),] %*% V),2,norm,type="2")
  order(omegas,decreasing = T)[1:WinningPoolSize]
  return(order(omegas,decreasing = T)[1:WinningPoolSize])
}


############# Factorization through coefficient matrix ######################
BPhyloFactor <- function(Data,tree,X,nfactors=1,formula=cbind(Successes,Failures)~y+z){
  getGLM <- function(spp,formula,DF){
    return(glm(formula,family=binomial,data=cbind(X,data.frame('Successes'=spp,'Failures'=1-spp))))
  }
  GLMs <- apply(Data,1,getGLM,formula=formula,DF=X)
  coefs <- sapply(GLMs,FUN=function(gg) coef(gg)[c('y','z')]/sqrt(diag(vcov(gg)))[c('y','z')])

  treeList <- list(tree)
  binList <- list(1:nrow(Data))
  Grps <- getPhyloGroups(tree)
  V <- sapply(Grps,ilrvec,nrow(Data))
  basis <- NULL
  
  for (pfs in 1:nfactors){
    if (pfs>1){
      treeList <- updateTreeList(treeList,binList,grp,tree,skip.check=T)
      binList <- updateBinList(binList,grp)
      Grps <- getNewGroups(tree,treeList,binList)
      V <- sapply(Grps,ilrvec,nrow(Data))
    }
    
    winner <- Bproj(V,coefs)
    grp <- Grps[[winner]]
    basis <- cbind(basis,ilrvec(grp,nrow(Data)))
  }
  return(basis)
}


#################### factorization through "phylo" factor ######################
## data tables will be used for quick referencing of species if mStableAgg=F
MatrixToDT <- function(l){
  Data <- l$Data
  X <- l$X
  DF <- data.table("Successes"=c(Data),
                   "Failures"=1-c(Data),
                   "Species"=rep(rownames(Data),times=ncol(Data)),
                   "Sample"=rep(colnames(Data),each=nrow(Data)),
                   key="Sample")
  X <- as.data.table(X)
  X[,Sample:=colnames(Data)]
  setkey(X,Sample)
  DF <- DF[X]
  setkey(DF,"Species")
  return(DF)
}

phyloFrame <- function(Data,grp,tree){
  factorFrame <- data.table('Species'=tree$tip.label,'phylo'='R')
  factorFrame[grp[[2]],'phylo'] <- 'S'
  return(Data[factorFrame])
}

mAgg <- function(Data,grp,tree,X){
  r <- length(grp[[1]])
  s <- length(grp[[2]])
  DF2 <- data.table('Sample'=rep(colnames(Data),times=2),
                    'Successes'=c(colSums(Data[grp[[1]],,drop=F]),colSums(Data[grp[[2]],,drop=F])),
                    'Failures'=rep(c(r,s),each=ncol(Data)),
                    'phylo'=factor(rep(c('R','S'),each=ncol(Data))),key='Sample')
  
  DF2 <- DF2[X]
  return(DF2)
}



glmPhyloFactor <- function(Data,tree,X=NULL,nfactors=1,mStableAgg=F){
  if (mStableAgg){
    X$Sample <- colnames(Data)
    X <- as.data.table(X)
    setkey(X,'Sample')
  }
  getDeviance <- function(fit){
    return(sum(anova(fit)[c('phylo:y','phylo:z'),'Deviance']))
  }
  getOmegas <- function(grp,Data,tree,X,mStableAgg,frmla=cbind(Successes,Failures)~phylo*(y+z)){
    if (mStableAgg){
      fit <- glm(formula=frmla,data=mAgg(Data,grp,tree,X),family=binomial)
    } else {
      fit <- glm(formula=frmla,data=phyloFrame(Data,grp,tree),family=binomial)
    }
    omega <- getDeviance(fit)
    return(omega)
  }
  
  treeList <- list(tree)
  binList <- list(1:(ape::Ntip(tree)))
  Grps <- getPhyloGroups(tree)
  basis <- NULL
  for (pfs in 1:nfactors){
    if (pfs>1){
      treeList <- updateTreeList(treeList,binList,grp,tree,skip.check=T)
      binList <- updateBinList(binList,grp)
      Grps <- getNewGroups(tree,treeList,binList)
    }
    
    omegas <- sapply(Grps,getOmegas,Data,tree,X,mStableAgg)
    winner <- which.max(omegas)
    grp <- Grps[[winner]]
    basis <- cbind(basis,ilrvec(grp,ape::Ntip(tree)))
  }
  
  return(basis)
}


############## mixed: uses B for quick computation, then picks the winning pool of edges for glmPhyloFactorization
mixPhyloFactor <- function(Data,tree,X,nfactors=1,WinningPoolSize=NULL,formula=cbind(Successes,Failures)~y+z){
  if (is.null(WinningPoolSize)){
    WinningPoolSize <- round(0.2*ape::Ntip(tree))
  }
  getGLM <- function(spp,formula,DF){
    return(glm(formula,family=binomial,data=cbind(X,data.frame('Successes'=spp,'Failures'=1-spp))))
  }
  getDeviance <- function(fit){
    return(sum(anova(fit)[c('phylo:y','phylo:z'),'Deviance']))
  }
  getOmegas <- function(grp,Data,frmla=cbind(Successes,Failures)~phylo*(y+z)){
    fit <- glm(formula=frmla,data=phyloFrame(Data,grp,tree),family=binomial)
    omega <- getDeviance(fit)
    return(omega)
  }
  GLMs <- apply(Data,1,getGLM,formula=formula,DF=X)
  D2 <- MatrixToDT(list('Data'=Data,'X'=X))
  coefs <- sapply(GLMs,FUN=function(gg) coef(gg)[c('y','z')]/sqrt(diag(vcov(gg)))[c('y','z')])

  treeList <- list(tree)
  binList <- list(1:nrow(Data))
  Grps <- getPhyloGroups(tree)
  V <- sapply(Grps,ilrvec,nrow(Data))
  basis <- NULL
  
  for (pfs in 1:nfactors){
    if (pfs>1){
      treeList <- updateTreeList(treeList,binList,grp,tree,skip.check=T)
      binList <- updateBinList(binList,grp)
      Grps <- getNewGroups(tree,treeList,binList)
      V <- sapply(Grps,ilrvec,nrow(Data))
    }
    
    tops <- Bproj(V,coefs,WinningPoolSize)
    
    omegas <- sapply(Grps[tops],getOmegas,D2)
    winner <- tops[which.max(omegas)]
    
    grp <- Grps[[winner]]
    basis <- cbind(basis,ilrvec(grp,nrow(Data)))
  }
  return(basis)
}


# Algorithm Comparison ----------------------------------------------------
v1 <- ilrvec(list(clade1,setdiff(1:50,clade1)),50) %>% sign
v2 <- ilrvec(list(clade2,setdiff(1:50,clade2)),50) %>% sign
Vref <- cbind(v1,v2)
## We can check for the correct edges by matching output basis vectors to v1 and v2.
## For instance, the vector below correctly partitions clade 2,
## after it has removed species 4 and 5 (presumably in a previous factor)
vv1 <- ilrvec(list(clade1,setdiff(1:50,clade1)),50)
vv2 <- ilrvec(list(clade2,setdiff(1:50,c(clade1,clade2))),50)

Vcheck <- function(V,Vref.=Vref){
  a <- apply(V,2,FUN=function(v) sum(v!=0))
  chk <- t(sign(V)) %*% Vref
  ix=0
  for (i in 1:ncol(V)){
    ix=ix+as.numeric(any(chk[i,]==a[i]))
  }
  return(ix)
}



reps=1000
correct <- c(0,0,0,0)
times <- c(0,0,0,0)

for (rr in 1:reps){
  l <- datasim()
  Data <- l$Data
  X <- l$X
  t1 <- Sys.time()
  Vbs <- BPhyloFactor(Data,tree,X,nfactors=2)
  t2 <- Sys.time()
  times[1] <- times[1]+difftime(t2,t1,units='secs')
  
  t1 <- Sys.time()
  VM <- glmPhyloFactor(Data,tree,X,mStableAgg=T,nfactors=2)
  t2 <- Sys.time()
  times[2] <- times[2]+difftime(t2,t1,units='secs')
  
  t1 <- Sys.time()
  Vtot <- glmPhyloFactor(MatrixToDT(l),tree,nfactors=2)
  t2 <- Sys.time()
  times[3] <- times[3]+difftime(t2,t1,units='secs')
  
  t1 <- Sys.time()
  Vmix <- mixPhyloFactor(Data,tree,X,nfactors=2)
  t2 <- Sys.time()
  times[4] <- times[4]+difftime(t2,t1,units='secs')
  
  correct <- correct + c(Vcheck(Vbs),
                         Vcheck(VM),
                         Vcheck(Vtot),
                         Vcheck(Vmix))
}

save(list=ls(),file='Algorithm_Comparison_workspace')



# Timing algorithms -------------------------------------------------------
set.seed(2)
mset <- c(50,100,150,200,250,300)
nfactors=1:3
algorithms <- c('B','mStable','phylo','mix')
scaling.times <- expand.grid(mset,nfactors,algorithms)
names(scaling.times) <- c('m','nfactors','algorithm')
scaling.times$time <- numeric(nrow(scaling.times))

for (m in mset){
  for (nf in nfactors){
    tree <- rtree(m)
    l <- datasim(n=40,m=m,tr=tree,c1=1,c2=2,cnst=rnorm(m))
    Data <- l$Data
    X <- l$X
    t1 <- Sys.time()
    invisible(BPhyloFactor(Data,tree,X,nfactors=nf))
    t2 <- Sys.time()
    scaling.times$time[scaling.times$m==m & scaling.times$algorithm=='B' & scaling.times$nfactors==nf] <- difftime(t2,t1,units='secs')
    
    t1 <- Sys.time()
    invisible(glmPhyloFactor(Data,tree,X,mStableAgg=T,nfactors=nf))
    t2 <- Sys.time()
    scaling.times$time[scaling.times$m==m & scaling.times$algorithm=='mStable' & scaling.times$nfactors==nf] <-difftime(t2,t1,units='secs')
    
    t1 <- Sys.time()
    invisible(glmPhyloFactor(MatrixToDT(l),tree,nfactors=nf))
    t2 <- Sys.time()
    scaling.times$time[scaling.times$m==m & scaling.times$algorithm=='phylo' & scaling.times$nfactors==nf] <-difftime(t2,t1,units='secs')
    
    t1 <- Sys.time()
    invisible(mixPhyloFactor(Data,tree,X,nfactors=nf))
    t2 <- Sys.time()
    scaling.times$time[scaling.times$m==m & scaling.times$algorithm=='mix' & scaling.times$nfactors==nf] <-difftime(t2,t1,units='secs')
  }
}



ggplot(scaling.times,aes(x=m,y=time,by=algorithm,color=algorithm))+
  geom_point(cex=5)+
  geom_line(lwd=2)+
  facet_grid(.~nfactors)+
  scale_y_continuous(trans='log',breaks = c(0.1,1,10,100,200,400))+
  scale_x_continuous(trans='log',breaks=mset)+
  theme_minimal()+
  theme(text = element_text(size=20),
        axis.text.x = element_text(angle=45, hjust=1)) 
# ggsave('gpf_algorithm_scaling.png',height=4,width=12)

save(list=ls(),file='Algorithm_Comparison_workspace')

getScalingCoef <- function(y,x) coef(glm(log(y)~log(x)))['log(x)']

scaling.times$coef <- numeric(nrow(scaling.times))
for (alg in algorithms){
  L <- scaling.times$algorithm==alg & scaling.times$nfactors==1
  y <- scaling.times$time[L]
  x <- scaling.times$m[L]
  scaling.times$coef[scaling.times$algorithm==alg] <- getScalingCoef(y,x)
}

getFactorCoef <- function(alg){
  L1 <- scaling.times$m==300 & scaling.times$algorithm==alg & scaling.times$nfactors==3
  L2 <- scaling.times$m==300 & scaling.times$algorithm==alg & scaling.times$nfactors==1
  return(scaling.times$time[L1]/scaling.times$time[L2])
}

scaling.times$factorScaling <- numeric(nrow(scaling.times))
for (alg in algorithms){
  scaling.times$factorScaling[scaling.times$algorithm==alg] <- getFactorCoef(alg)
}
# Algorith comparison constant mean ---------------------------------------
### Now we compare algorithms when the mean is constant.

set.seed(1)
m=50 ## number of species
n=40 ## number of samples
tree <- rtree(m)

reps2=3000
correct2 <- c(0,0,0,0)

for (rr in 1:reps2){
  l <- datasim(cnst=rep(0,m))
  Data <- l$Data
  X <- l$X
  
  Vbs <- BPhyloFactor(Data,tree,X,nfactors=2)
  VM <- glmPhyloFactor(Data,tree,X,mStableAgg=T,nfactors=2)
  Vtot <- glmPhyloFactor(MatrixToDT(l),tree,nfactors=2)
  Vmix <- mixPhyloFactor(Data,tree,X,nfactors=2)
  
  correct2 <- correct2 + c(Vcheck(Vbs),
                         Vcheck(VM),
                         Vcheck(Vtot),
                         Vcheck(Vmix))
}

save(list=ls(),file='Algorithm_Comparison_workspace')
