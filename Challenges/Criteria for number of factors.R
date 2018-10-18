#### Stopping Criteria

## In this script, we will investigate the following stopping criteria:
# KS-test of P-values from regression,
# Horn-style factorization of null dataset (%ex var)
# row & column re-sampling of null datsaet (%ex var)

## To do this, we will consider a 32-species tree and randomly assign 
## effects to 2, 4, 8, or 16 nodes (excluding root). 

## Then, we will consider two figures: 

# (1) The % Ex var under each simulation (excluding KS-test)
       #We'll use this to motivate a stopping criteria for each of Horn & row/col resampling
# (2) The true vs. estimated # of factors

# These results may be sensitive to tree-topology, but we will start off considering only one tree.
library(ape)
library(phytools)
library(phylofactor)
library(parallel)
setwd('Challenges')
set.seed(1)
nspecies=32
tree <- rtree(nspecies)
PAR <- par('mar')

clr <- function(A) apply(A,MARGIN=2,FUN=function(A)log(A)-mean(log(A)))


effect.sim <- function(n.effects=2,m=10,X=NULL,tree,effects=3:18-0.1){
  ### Initialize dataset
  n <- Ntip(tree)
  Data <- rlnorm(n*m) %>% matrix(.,nrow=n)
  rownames(Data) <- tree$tip.label
  if (is.null(X)){X=rnorm(m)}
  
  ### Draw clades
  nodes <- sample(c(1:n,(n+3):(2*n-1)),n.effects,replace=F) #skipping n+1 and n+2 avoids root and its lower descendant
  tips <- nodes[nodes<=n]
  nodes <- setdiff(nodes,tips)
  if (length(nodes)>0){
    clades <- lapply(as.list(nodes),FUN=function(node,tree) phangorn::Descendants(tree,node,'tips')[[1]],tree=tree)
  } else {clades <- NULL}
  clades <- c(tips,clades)
  
  ### Draw effects
  effects <- sample(effects,n.effects,replace=F)*sign(rnorm(n.effects))

  ### assign effects
  for (nn in 1:n.effects){
    cl <- clades[[nn]]
    for (kk in 1:length(cl)){
      Data[cl[kk],] <- Data[cl[kk],]*exp(effects[nn]*X)
    }
  }
  
  return(list('Data'=Data,'X'=X))
}

# Demo:
# Y <- effect.sim(tree=tree)
# phylo.heatmap(tree,clr(Y$Data[,order(Y$X)]))

#obtain replicate Explained-Variance matrices
repEV <- function(n.effects,reps=100,nspecies=32,type='pf',method='Shuffle'){
  
  if (! type %in% c('pf','phyca','null')){
    stop('input type must be either pf or phyca')
  }
  Ex.var <- matrix(NA,nrow=reps,ncol=nspecies-1)
  for (i in 1:reps){
    Y <- effect.sim(n.effects = n.effects,tree=tree)
    if (type=='pf'){
      pf <- PhyloFactor(Y$Data,tree,Y$X)
      Ex.var[i,] <- pf$ExplainedVar
    } else if (type=='phyca') {
      phyca <- PhyCA(Y$Data,tree,ncomponents=nspecies-1)
      Ex.var[i,] <- phyca$PercentVariance
    } else {
      pf <- list('Data'=Y$Data,'tree'=tree,'X'=Y$X)
      pf <- pf.nullsim(pf,1,method=method,nfactors=nspecies-1)
      Ex.var[i,] <- pf[[1]]
    }
  }
  return(Ex.var)
}

REPS=300

cl <- phyloFcluster(ncores=4)
clusterExport(cl,varlist=c('repEV','effect.sim','tree'))
Neffects <- c(2,4,8,16)
EV <- parLapply(cl,as.list(Neffects),fun=repEV,reps=REPS)
names(EV) <- Neffects
stopCluster(cl)
rm('cl')




tiff('ExVar_curves_Neffects_2-8.tiff',height=1000,width=1000)
  par(mfrow=c(2,2))
  for (n in Neffects){
    Ex.var <- EV[toString(n)][[1]]
    for (i in 1:nrow(Ex.var)){
      if (i==1){
        plot(1:(nspecies-1),Ex.var[i,],type='l',log='y',ylim=c(1e-7,1),main=paste(n,'affected clades'),col='grey',cex.main=2,cex.lab=2,cex.axis=2)
        lines(1:(nspecies-1),Ex.var[i,],col='grey')
      } else {
        lines(Ex.var[i,],col='grey')
      }
    }
    ml <- colMeans(Ex.var)
    lines(1:(nspecies-1),ml,lwd=4)
    lines(c(n,n),c(1e-10,2),lty=2,lwd=4,col='red')
  }
dev.off()

save(list=ls(),file='Criteria_for_number_of_factors_workspace')


## two kinds of null simulations:

## "Gaussian" null simulations - null data are i.i.d. standard log-normal, so we'll use the same 100 reps for all Neffects 
Y <- effect.sim(tree=tree)
null.Gaussian <- pf.nullsim(PF=list('Data'=Y$Data,'tree'=tree,'X'=Y$X),reps=REPS)
NG <- unlist(null.Gaussian) %>% matrix(.,nrow=REPS,byrow=T)

## "Shuffle" null simulations - shuffles rows and columns of dataset
cl <- phyloFcluster(ncores=4)
clusterExport(cl,varlist=c('repEV','effect.sim','tree','nfactors'))
Neffects <- c(2,4,8,16)
null.Shuffle <- parLapply(cl,as.list(Neffects),fun=repEV,type='null',reps=REPS)
names(null.Shuffle) <- Neffects
stopCluster(cl)
rm('cl')

save(list=ls(),file='Criteria_for_number_of_factors_workspace')



# 
# tiff('ExVar_curves_with_nullSims.tiff',height=1000,width=1000)
# par(mfrow=c(2,2))
# for (n in Neffects){
#   Ex.var <- EV[toString(n)][[1]] 
#   NS <- null.Shuffle[toString(n)][[1]]
#   ml.ns <- colMeans(NS)
#   ml.ng <- colMeans(NG)
#   
#   ### Initialize plot
#   plot(1:(nspecies-1),Ex.var[1,],type='l',log='y',ylim=c(1e-7,1),main=paste(n,'affected clades'),col='grey',xlab='factor',ylab='% Variance Explained',cex.main=2,cex.lab=2,cex.axis=2)
#   
#   
#   ### Plot Gaussian null simulations
#   for (i in 1:nrow(NS)){
#     lines(1:(nspecies-1),null.Gaussian[[i]],col=rgb(0,1,0,alpha=0.1))
#   }
#   lines(1:(nspecies-1),ml.ng,col='green',lwd=4)
#   
#   ### Plot Shuffle null simulations
#   for (i in 1:nrow(NS)){
#     lines(1:(nspecies-1),NS[i,],col=rgb(0,0,1,alpha=0.1)) 
#   }
#   lines(1:(nspecies-1),ml.ng,col='green',lwd=4,pch=1,type='o',cex=2)
#   lines(1:(nspecies-1),ml.ns,col='blue',lwd=4,pch=4,type= 'o',cex=2)
#   
#   
#   ml <- colMeans(Ex.var)
#   lines(1:(nspecies-1),ml,lwd=4)
#   lines(c(n,n),c(1e-10,2),lty=3,lwd=4,col='red')
#   if (n==2){
#     legend('topright',c('Phylofactor Mean','null: log-normal','null: shuffle','true # factors'),col=c('black','green','blue','red'),lty=c(1,1,1,3),pch=c(NA,1,4,NA),lwd=c(4,4,4,4),cex=2)
#   }
# }
# dev.off()


######################### Stopping Criteria ##################

####### null-LN

nfactsML=function(ev,ml.ng)  apply(ev,MARGIN=1,FUN=function(ev,ml.ng) min(which((ev-ml.ng)<0)),ml.ng=ml.ng)
NF <- lapply(EV,nfactsML,ml.ng=ml.ng) %>% as.data.frame
names(NF) <- Neffects




####### KS
## For comparison, we need KS stopping points for the same.
KS.stop <- function(n.effects,reps,tree,Pval=0.01,alternative='two.sided'){
  Nfactors <- rep(NA,reps)
  for (rr in 1:reps){
    Y <- effect.sim(n.effects = n.effects,tree=tree)
    pf <- PhyloFactor(Y$Data,tree,Y$X,stop.early=F,KS.Pthreshold = Pval,alternative=alternative)
    Nfactors[rr] <- pf$nfactors
  }
  return(Nfactors)
}

cl <- phyloFcluster(ncores=4)
clusterExport(cl,varlist=c('tree','effect.sim','KS.stop'))
KS0.01 <- parLapply(cl,as.list(Neffects),KS.stop,reps=REPS,tree=tree)
stopCluster(cl)
rm('cl')

names(KS0.01) <- Neffects

cl <- phyloFcluster(ncores=4)
clusterExport(cl,varlist=c('tree','effect.sim','KS.stop'))
KSgreater <- parLapply(cl,as.list(Neffects),KS.stop,reps=REPS,tree=tree,alternative='greater')
stopCluster(cl)
rm('cl')

names(KSgreater) <- Neffects

save(list=ls(),file='Criteria_for_number_of_factors_workspace')



####### THE figure:


tiff('KS vs LN stopping criteria.tiff',height=1000,width=2000)
  layout(matrix(1:8,nrow=2,byrow=F))
  for (n in Neffects){
    Ex.var <- EV[toString(n)][[1]] 
    ml.ng <- colMeans(NG)
    
    ### Initialize plot
    ml <- colMeans(Ex.var)
    plot(1:(nspecies-1),ml,type='l',log='y',ylim=c(1e-7,1),main=paste(n,'affected clades'),col='grey',xlab='factor',ylab='% Variance Explained',cex.main=3,cex.lab=3,cex.axis=2)
    
    
    ### Plot Gaussian null simulations
    for (i in 1:nrow(NS)){
      lines(1:(nspecies-1),null.Gaussian[[i]],col=rgb(0,1,0,alpha=0.05))
    }
    lines(1:(nspecies-1),ml.ng,col='green',lwd=4)
    lines(1:(nspecies-1),ml.ng,col='green',lwd=4,pch=1,type='o',cex=2)

    
    
    lines(1:(nspecies-1),ml,lwd=4)
    lines(c(n,n),c(1e-10,2),lty=3,lwd=4,col='red')
    if (n==2){
      legend('topright',c('Phylofactor Mean','Log-normal null','true # factors'),col=c('black','green','red'),lty=c(1,1,3),pch=c(NA,1,NA),lwd=c(4,4,4),cex=3)
    }
    
    DF <- data.frame('LN'=NF[,toString(n)]-1,'KS'=KSgreater[[toString(n)]]-1)
    boxplot(DF,ylim=c(n/2,2*n),cex=3,lwd=2,cex.lab=3,cex.axis=3,ylab='# factors')
    hrn <- c(max(DF$LN),var(DF$LN),sum(DF$LN>n)/300) %>% sapply(.,signif,digits=3) %>% mapply(FUN=function(a,b) paste(b,a,sep=''),b=c('max=','var=','OF Rate='))
    ks <- c(max(DF$KS),var(DF$KS),sum(DF$KS>n)/300) %>% sapply(.,signif,digits=3) %>% mapply(.,FUN=function(a,b) paste(b,a,sep=''),b=c('max=','var=','OF Rate='))
    lines(c(-10,10),c(n,n),col='red',lty=3,lwd=4)
    
    
    legend('topleft',hrn,cex=3)
    legend('topright',ks,cex=3)
  }
dev.off()


library(ggplot2)
library(ggpubr)
cols <- c('green',viridis(1))


n=2
DF <- data.table('factor'=c(NF[,toString(n)]-1,KSgreater[[toString(n)]]-1),'cutoff'=rep(c('Horn','KS'),each=nrow(NF)))
x <- sort(unique(DF$factor))
cnt <- max(DF[,table(factor),by=cutoff]$V1)
g2 <- ggplot(DF,aes(factor,fill=cutoff,color=cutoff))+
  geom_histogram(aes(y=2*..count../(sum(..count..))),position='identity',alpha=0.3,binwidth=1,center=0)+
  geom_vline(xintercept = n+.5,col='red',lwd=3,lty=2)+
  scale_fill_manual(values=cols,labels=c("Horn","KS"))+
  scale_y_continuous('')+
  theme_minimal()+
  theme(legend.position = c(0.75,.8),
        axis.text.y = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.title = element_text(size=0),
        legend.key.size = unit(x=.4,units = 'in'),
        text = element_text(size=20))
n=4
DF <- data.table('factor'=c(NF[,toString(n)]-1,KSgreater[[toString(n)]]-1),'cutoff'=rep(c('Horn','KS'),each=nrow(NF)))
x <- sort(unique(DF$factor))
cnt <- max(DF[,table(factor),by=cutoff]$V1)
g4 <- ggplot(DF,aes(factor,fill=cutoff,color=cutoff))+
  geom_histogram(aes(y=2*..count../sum(..count..)),position='identity',alpha=0.3,binwidth=1,center=0)+
  geom_vline(xintercept = n+.5,col='red',lwd=3,lty=2)+
  scale_fill_manual(values=cols,labels=c("Horn","KS"))+
  scale_y_continuous('')+
  theme_minimal()+
  theme(legend.position = c(0.75,.8),
        axis.text.y = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.title = element_text(size=0),
        legend.key.size = unit(x=.4,units = 'in'),
        text = element_text(size=20))
n=8
DF <- data.table('factor'=c(NF[,toString(n)]-1,KSgreater[[toString(n)]]-1),'cutoff'=rep(c('Horn','KS'),each=nrow(NF)))
x <- sort(unique(DF$factor))
cnt <- max(DF[,table(factor),by=cutoff]$V1)
g8 <- ggplot(DF,aes(factor,fill=cutoff,color=cutoff))+
  geom_histogram(aes(y=2*..count../sum(..count..)),position='identity',alpha=0.3,binwidth=1,center=0)+
  geom_vline(xintercept = n+.5,col='red',lwd=3,lty=2)+
  scale_fill_manual(values=cols,labels=c("Horn","KS"))+
  scale_x_continuous(breaks=c(3,6,9,12))+
  scale_y_continuous('')+
  theme_minimal()+
  theme(legend.position = c(0.2,.8),
        axis.text.y = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.title = element_text(size=0),
        legend.key.size = unit(x=.4,units = 'in'),
        text = element_text(size=20))
n=16
DF <- data.table('factor'=c(NF[,toString(n)]-1,KSgreater[[toString(n)]]-1),'cutoff'=rep(c('Horn','KS'),each=nrow(NF)))
x <- sort(unique(DF$factor))
cnt <- max(DF[,table(factor),by=cutoff]$V1)
g16 <- ggplot(DF,aes(factor,fill=cutoff,color=cutoff))+
  geom_histogram(aes(y=2*..count../sum(..count..)),position='identity',alpha=0.3,binwidth=1,center=0)+
  geom_vline(xintercept = n+.5,col='red',lwd=3,lty=2)+
  scale_fill_manual(values=cols,labels=c("Horn","KS"))+
  scale_y_continuous('')+
  theme_minimal()+
  theme(legend.position = c(0.2,.8),
        axis.text.y = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.title = element_text(size=0),
        legend.key.size = unit(x=.4,units = 'in'),
        text = element_text(size=20))

ggarrange(g2,g4,g8,g16,ncol=4)
ggsave('KS_LN_compare.png',width=18,height=4.5)

