####### Phylofactorization as unbiased graph-partitioning algorithm
library(phylofactor)
library(parallel)
library(RMTstat) ## for Marchenko-Pastur distribution
library(EMT)
library(phytools)
library(ggplot2)
library(reshape2)
setwd('Challenges')



# Functions ---------------------------------------------------------------

### grab's an edge's distance to root
distance.to.root <- function(edge){
  nd <- tree$edge[edge,1]
  np <- nodepath(tree,nd,Ntip(tree)+1)
  eln <- 0
  for (i in 1:(length(np)-1)){
    ed <- apply(tree$edge,1,FUN=function(a,b) all(a%in%b),b=np[c(i,i+1)])
    eln <- eln+tree$edge.length[which(ed)]
  }
  return(eln)
}

### simulate's null data
null.data <- function(D=10,n=8,tree=tree,meanlog=0,sdlog=1,Xtype='real',rvar='lnorm'){
  if (rvar=='lnorm'){
    Data <- matrix(rlnorm(n=D*n,meanlog=meanlog,sdlog=sdlog),nrow=D)
  } else if (rvar=='fastBM'){
    Data <- exp(sapply(1:n,FUN=function(n,tree) fastBM(tree = tree),tree=tree))
  } else {
    Data <- matrix(rnbinom(n=D*n,size=1,prob = 5e-4),nrow=D)
  }
  rownames(Data) <- tree$tip.label
  if (Xtype=='real'){
    X <- rnorm(n)
  } else {
    X <- rep(c(1,2),each=n/2)
  }
  return(list('Data'=Data,'X'=X))
}

### computes one phylofactor & gets factored edge
pf.edges <- function(reps,tree=ape::rtree(20),output.gdiff=F,n=8,choice='F',method='glm',meanlog=0,sdlog=1,Xtype='real',rvar='lnorm'){
  edgs <- vector(mode='list',length=reps)
  D <- length(tree$tip.label)
  for (rr in 1:reps){
    Dum <- null.data(D=D,n=n,tree=tree,meanlog=meanlog,sdlog=sdlog,Xtype=Xtype,rvar=rvar)
    Dum$Data[Dum$Data==0]=0.65
    PF <- PhyloFactor(Data=Dum$Data,tree=tree,X=Dum$X,choice=choice,method=method,nfactors=1)
    edgs[[rr]] <- getFactoredEdges(v=PF$basis[,1,drop=F],tree=tree)
  }
  return(edgs)
}

### Fully phylofactors a null dataset for method='max.var'
### returns edge for first factor + explained var for each factor
pf.maxvar <- function(reps,tree=ape::rtree(10),n=8,choice='var',meanlog=0,sdlog=1,Xtype='real',rvar='lnorm'){
  edgs <- vector(mode='list',length=reps)
  exvar <- edgs
  D <- length(tree$tip.label)
  
  for (rr in 1:reps){
    Dum <- null.data(D=D,n=n,tree=tree,meanlog=meanlog,sdlog=sdlog,Xtype=Xtype,rvar=rvar)
    PF <- PhyloFactor(Data=Dum$Data,tree=tree,X=Dum$X,choice=choice,method='max.var',nfactors=(D-1))
    edgs[[rr]] <- getFactoredEdges(v=PF$basis[,1,drop=F],tree=tree)
    exvar[[rr]] <- PF$ExplainedVar*PF$total.variance
  }
  
  return(list('edges'=edgs,'ExplainedVar'=exvar))
}


#### random objective function phylofactorization
pf.rand <- function(reps,tree){
  D <- Ntip(tree)
  Data <- matrix(rlnorm(8*D),nrow=D)
  rownames(Data) <- tree$tip.label
  X <- rnorm(8)
  
  ro <- function(y,X,PF.output=FALSE){
    output <- NULL
    if(PF.output){
      return(NA)
    } else {
      output$objective <- rlnorm(1)
      output$stopStatistics <- 0
      return(output)
    }
  }
  
  edgs <- vector(mode='list',length=reps)
  for (rr in 1:reps){
    PF <- PhyloFactor(Data,tree,X,choice.fcn = ro,nfactors=1)
    edgs[[rr]] <- getFactoredEdges(v=PF$basis[,1,drop=F],tree=tree)
  }
  return(edgs)
}




# Random phylofactorization ---------------------------------------------------------


set.seed(1)
tree=rtree(40)
tree$edge.length=rep(0.1,Nedge(tree))
ncores=7
cl <- phyloFcluster(ncores)
clusterExport(cl,varlist='pf.rand')
reps=16000
reps <- as.list(rep(reps,ncores))
edgs <- parLapply(cl,reps,pf.rand,tree=tree)

set.seed(2)

edgs <- unlist(edgs)
edgR <- edgs[!edgs==1]   #remove one root edge
tb <- table(edgR)

### visualize this by weighting edges
els <- sapply(as.numeric(names(tb)),distance.to.root)
root.edge <- which(sapply(els,length)==0)
els[root.edge]=0
w <- tb/min(tb)
wR <- c(w[root.edge],w)

plot.phylo(tree,edge.width = 2*wR^2)
chisq.test(tb)   #p=0.4145






# Objective-functions for phylofactorization------------------------------------------
D <- Ntip(tree)
E <- 2*D-3

root.edges <- which(tree$edge[,1]==(D+1))
tips <- which(tree$edge[,2] %in% 1:D)
clusterExport(cl,varlist=c('null.data','pf.edges','pf.maxvar','tree','Rgpf'))



# Var ---------------------------------------------------------------------

edgs <- parLapply(cl,X=reps,fun=pf.edges,tree=tree,choice='var')

edgs <- unlist(edgs)
edgV <- edgs[!edgs==1]   #remove one root edge
tb <- table(edgV)

### visualize this by weighting edges
els <- sapply(as.numeric(names(tb)),distance.to.root)
root.edge <- which(sapply(els,length)==0)
els[root.edge]=0
w <- tb/min(tb)
wV <- c(w[root.edge],w)


# Negative Binomial Counts ------------------------------------------------

#### rnbinom counts
edgs <- parLapply(cl,X=reps,fun=pf.edges,tree=tree,rvar='nbinom')

edgs <- unlist(edgs)
edgNB <- edgs[!edgs==1]   #remove one root edge
tb <- table(edgNB)

#visualize 
els <- sapply(as.numeric(names(tb)),distance.to.root)
root.edge <- which(sapply(els,length)==0)
els[root.edge]=0
w <- tb/min(tb)
wNB <- c(w[root.edge],w)


### rnbinom choice='var
edgs <- parLapply(cl,X=reps,fun=pf.edges,tree=tree,rvar='nbinom',choice='var')

edgs <- unlist(edgs)
edgNBV <- edgs[!edgs==1]   #remove one root edge
tb <- table(edgNBV)

#visualize 
els <- sapply(as.numeric(names(tb)),distance.to.root)
root.edge <- which(sapply(els,length)==0)
els[root.edge]=0
w <- tb/min(tb)
wNBV <- c(w[root.edge],w)



# F -----------------------------------------------------------------------

###### choice='F'
edgs <- parLapply(cl,X=reps,fun=pf.edges,tree=tree,choice='F')

edgs <- unlist(edgs)
edgF <- edgs[!edgs==1]   #remove one root edge
tb <- table(edgF)
els <- sapply(as.numeric(names(tb)),distance.to.root)
root.edge <- which(sapply(els,length)==0)
els[root.edge]=0
### visualize this by weighting edges
w <- tb/min(tb)
wF <- c(w[root.edge],w)






# max.var -----------------------------------------------------------------
edgs <- parLapply(cl,X=reps,fun=pf.edges,tree=tree,choice='var',method='max.var')

edgs <- unlist(edgs)
edgMV <- edgs[!edgs==1]   #remove one root edge
tbMV <- table(edgMV)
els <- sapply(as.numeric(names(tb)),distance.to.root)
root.edge <- which(sapply(els,length)==0)
els[root.edge]=0
### visualize this by weighting edges
w <- tb/min(tb)
wMV <- c(w[root.edge],w)

stopCluster(cl)
rm('cl')

# save(list=ls(),file='graph_partitioning_analysis_big_workspace_2')



# fastBM ------------------------------------------------------------------
cl <- phyloFcluster(ncores)
clusterExport(cl,varlist=c('pf.edges','null.data'))
clusterEvalQ(cl,library(phytools))
edgs <- parLapply(cl,reps,pf.edges,tree=tree,rvar='fastBM')
stopCluster(cl)
rm('cl')


edgs <- unlist(edgs)
edgBM <- edgs[!edgs==1]   #remove one root edge
tb <- table(edgBM)

# save(list=ls(),file='graph_partitioning_analysis_big_workspace_2')



# Plotting evenness vs P(Edge) --------------------------------------------
# load('graph_partitioning_analysis_big_workspace_2')


clo <- function(x) x/sum(x)
tbR <- table(edgR) %>% clo
tbMV <- table(edgMV) %>% clo
tbF <- table(edgF) %>% clo
tbNB <- table(edgNB) %>% clo
tbBM <- table(edgBM) %>% clo
TB <- tbMV+tbF+tbNB %>% clo
els <- sapply(as.numeric(names(tbR)),distance.to.root)
edg.evn <- function(edge,tree.=tree){
  r <- length(phangorn::Descendants(tree,tree$edge[edge,2],type='tips')[[1]])
  return(abs(Ntip(tree)-2*r))
}
gdiff <- as.factor(sapply(as.numeric(names(tbR)),edg.evn))

DF <- data.frame('edgeProb'=c(tbR,tbMV,tbF,tbNB,tbBM),'ObjFcn'=rep(c('Random','PVA','F-stat','F-stat NBinom','GraphBM'),each=length(tbR)),'Unevenness'=rep(gdiff,5))

DF <- DF[!DF$ObjFcn=='Random',]
sings <- which(table(DF$Unevenness)==4) %>% names %>% as.numeric
DF <- DF[!DF$Unevenness %in% sings,]
df <-data.frame('x'=c(-Inf,Inf),'y'=c(1/E,1/E))

DF$ObjFcn <- relevel(DF$ObjFcn,ref = 'GraphBM')
DF$ObjFcn <- factor(DF$ObjFcn,levels=rev(levels(DF$ObjFcn)))


par(mfrow=c(1,1))
DF$Unevenness <- factor(DF$Unevenness)
ggplot(DF,aes(x=Unevenness,y=edgeProb))+
  geom_boxplot(aes(fill=ObjFcn))+
  geom_line(data=df,aes(x=x,y=y))+
  ggtitle('Bias of ILR Graph-partitioning')+
  scale_y_continuous(name='Probability of Drawing Edge')+
  scale_x_discrete(name='Unevenness, |r-s|')+
  theme(title = element_text(size=20))

names(df) <- c('ex','why')
DF$Unevenness <- as.numeric(sapply(DF$Unevenness,toString))
# DF$Unevenness <- 
ggplot(DF,aes(x=Unevenness,y=edgeProb,fill=ObjFcn,color=ObjFcn)) + 
  geom_point(data=DF,size=3)+
  scale_y_continuous(name='Probability of Drawing Edge')+
  scale_x_discrete(name='Partition Unevenness, |r-s|')+
  theme(title = element_text(size=20))+
  facet_grid(.~ObjFcn)+
  geom_smooth(data=DF,method='glm',formula = y~exp(x))

ggsave(filename = '~/PhyloFactor/generalized_phylofactorization/Figures/PEdge_vs_Unevenness_by_Obj_Fcn.pdf',height = 4,width=16)
  
unique(DF$ObjFcn)
ddf <- data.frame('PVA' = DF$edgeProb[DF$ObjFcn=='PVA'],
                  'F-stat' = DF$edgeProb[DF$ObjFcn=='F-stat'],
                  'F-stat NBinom' = DF$edgeProb[DF$ObjFcn=='F-stat NBinom'],
                  'GraphBM' = DF$edgeProb[DF$ObjFcn=='GraphBM'])




# save(list=ls(),file='graph_partitioning_analysis_big_workspace_3')




# Voronoi cell estimation -------------------------------------------------
# load('graph_partitioning_analysis_big_workspace_3')
Grps <- getPhyloGroups(tree)
V <- sapply(Grps,ilrvec,D)
edgs <- apply(V,2,getFactoredEdges,tree=tree)
root.grps <- which(sapply(edgs,length)>1)

V <- V[,2:ncol(V)] %>% t
edgs <- edgs[2:length(edgs)]
Grps <- Grps[2:length(Grps)]
edgs[sapply(edgs,length)>1] <- root.edges[2]
edgs <- unlist(edgs)
nrow(V)==length(edgs)

### to simulate Voronoi cell estimation, we will simulate gaussian random vectors, y, and find
### which row of V, v, maximizes the projection of y onto v

Vmax <- function(reps,V,edgs){
  A <- matrix(rnorm(reps*ncol(V)),nrow=ncol(V))
  mx <- (V %*% A) %>% apply(.,MARGIN=2,FUN=function(x) which.max(abs(x)))
  tb <- table(edgs[mx])
  return(tb)
}

repsV <- as.list(rep(1e6,7))
cl <- makeCluster(7)
clusterExport(cl,varlist='Vmax')
clusterEvalQ(cl,library(magrittr))
TB <- parLapply(cl,repsV,Vmax,V=V,edgs=edgs)
stopCluster(cl)
rm('cl')
gc()

tbVoronoi <- TB[[1]]
for (i in 2:length(repsV)){
  tbVoronoi <- tbVoronoi + TB[[i]]
}
tbVoronoi <- clo(tbVoronoi)


ddf <- data.frame('PVA'=c(tbMV),'F-stat'=c(tbF),'NB F-stat'=c(tbNB),
                  'PVA-GraphBM'=c(tbBM),'Voronoi'=c(tbVoronoi))




tiff('~/PhyloFactor/generalized_phylofactorization/Figures/EdgeBiase_across_ObjFcns.tiff',height=800,width=800)
  plot(ddf,pch=16)
dev.off()



# save(list=ls(),file='graph_partitioning_analysis_big_workspace_4')

# Fisher Factor bias ------------------------------------------------------


fishtest <- function(reps,tree,Grps,edgs,p=0.5){
  Ftest <- function(grps,tree,Z){
    s1 <- Z[grps[[1]]]
    s2 <- Z[grps[[2]]]
    n1 <- sum(s1)
    n2 <- sum(s2)
    ft <- fisher.test(matrix(c(n1,length(s1)-n1,n2,length(s2)-n2),ncol=2),alternative='two.sided')
    return(log(ft$p.value))
  }
  TestFunction <- function(Z,tree,Grps){
    return(sapply(Grps,Ftest,tree=tree,Z=Z))
  }
  n <- Ntip(tree)
  E <- numeric(reps)
  
  for (i in 1:reps){
    Z <- rbinom(n,1,p)
    E[i] <- TestFunction(Z,tree,Grps) %>% which.min %>% edgs[.]
  }
  return(table(E))
}

repsF <- as.list(rep(1e4,7))
p=0.5

cl <- phyloFcluster(7)
clusterExport(cl,'fishtest')

TB.f <- parLapply(cl,repsF,fishtest,tree=tree,Grps=Grps,edgs=edgs)

tbFish <- tbVoronoi
tbFish[1:length(tbFish)] <- 0  #this gets the names right - not all edges were selected

for (i in 1:length(TB.f)){
  tt <- TB.f[[i]]
  tt <- tt[match(names(tbFish),names(tt))]
  names(tt) <- names(tbFish)
  tt[is.na(tt)]=0
  tbFish <- tbFish+tt
}



DF <- rbind(DF,data.frame('edgeProb'=c(tbFish)/sum(tbFish),'ObjFcn'='FisherTest','Unevenness'=gdiff))


DFF <- DF
DFF <- rbind(data.frame('edgeProb'=c(tbR)/sum(tbR),'ObjFcn'='Random','Unevenness'=gdiff),DFF)
DFF$edgeProb[DFF$edgeProb==0] <- min(DFF$edgeProb[DFF$edgeProb>0])

DFF$Unevenness <- as.numeric(DFF$Unevenness)
ggplot(DFF,aes(x=Unevenness,y=edgeProb,fill=ObjFcn,color=ObjFcn)) + 
  geom_point(data=DFF,size=3)+
  scale_y_continuous(name='Probability of Drawing Edge',trans='log')+
  scale_x_discrete(name='Partition Unevenness, |r-s|')+
  theme(title = element_text(size=20))+
  facet_grid(.~ObjFcn)+
  geom_smooth(data=DFF,method='glm',formula = y~exp(x))

ggsave(filename = '~/PhyloFactor/generalized_phylofactorization/Figures/PEdge_vs_Unevenness_by_Obj_Fcn_ALL.pdf',height = 3,width=16)


# save(list=ls(),file='graph_partitioning_analysis_big_workspace_final')






# gpf ---------------------------------------------------------------------

######### generalized phylofactorization for Presence/Absence data #########
Rgpf <- function(reps=1,tree,Samples=as.character(1:20),prob=0.5){
  
  ###### species, groups & edges
  spp <- tree$tip.label
  n <- length(Samples)
  D=length(spp)
  edgs <- numeric(reps)
  Grps <- phylofactor::getPhyloGroups(tree)
  V <- sapply(Grps,phylofactor::ilrvec,n=D)
  Edges <- apply(V,2,phylofactor::getFactoredEdges,tree=tree)
  root.edges <- which(sapply(Edges,length)==2)
  Edges[root.edges] <- root.edges[2]
  
  Grps <- Grps[2:length(Grps)]
  Edges <- unlist(Edges[2:length(Edges)])
  SampleFrame  <- expand.grid(Samples,c('R','S'))
  names(SampleFrame) <- c('Samples','G')
  SampleFrame <- as.data.table(SampleFrame)
  binom.size=1
  frmla=cbind(Successes,Failures)~zz*phylo
  expfamily='binomial'
  PartitioningVariables='zz'
  model.fcn=glm
  for (rr in 1:reps){
    PA <- lapply(spp,FUN=function(x,Samples,prob) Samples[as.logical(rbinom(n,binom.size,prob))],Samples,prob)
    M <- data.table('Species'=rep(spp,sapply(PA,length)),'Samples'=unlist(PA),'N'=1,stringsAsFactors = F) %>%
      phylo.frame.to.matrix
    X <- data.table('Samples'=Samples,'zz'=rnorm(n))
    ix <- which.max(sapply(Grps,getObjective,tree,M,X,binom.size,frmla,expfamily,model.fcn,PartitioningVariables,family=binomial))
    edgs[rr] <- Edges[ix]
  }
  return(edgs)
}


cl <- phyloFcluster(7)
parallel::clusterEvalQ(cl,library(data.table))
parallel::clusterExport(cl,varlist=c('Rgpf'))


edgGPF <- parallel::parSapply(cl,reps,FUN=Rgpf,tree=tree)

parallel::stopCluster(cl)
rm('cl')

tbGPF <- table(edgGPF)
tbGPF
chisq.test(tbGPF)


plot(c(table(edgV)),c(table(edgGPF)))


# DF <- rbind(DF,data.frame('edgeProb'=c(tbFish)/sum(tbFish),'ObjFcn'='FisherTest','Unevenness'=gdiff))
# 
# 
DFF <- DF
DFF <- rbind(data.frame('edgeProb'=c(tbR)/sum(tbR),'ObjFcn'='Random','Unevenness'=gdiff),DFF)
DFF$edgeProb[DFF$edgeProb==0] <- min(DFF$edgeProb[DFF$edgeProb>0])

DFF <- rbind(data.frame('edgeProb'=c(tbGPF)/sum(tbGPF),'ObjFcn'='GPF','Unevenness'=gdiff),DFF)
DFF <- DFF[DFF$ObjFcn!='F-stat NBinom',]
DFF$ObjFcn <- factor(DFF$ObjFcn,levels = c('Random','PVA','F-stat','GPF','GraphBM','FisherTest'))
DFF$Unevenness <- as.numeric(DFF$Unevenness)
ggplot(DFF,aes(x=Unevenness,y=edgeProb,fill=ObjFcn,color=ObjFcn)) + 
  geom_point(data=DFF,size=3)+
  scale_y_continuous(name='Probability of Drawing Edge',trans='log')+
  scale_x_discrete(name='Partition Unevenness, |r-s|')+
  theme(title = element_text(size=20))+
  facet_grid(.~ObjFcn)+
  geom_smooth()
ggsave(filename = 'PEdge_vs_Unevenness_by_Obj_Fcn_NEW.pdf',height = 3,width=16)
ggsave(filename = 'PEdge_vs_Unevenness_by_Obj_Fcn_NEW.tiff',height = 3,width=16)


ddf <- data.frame('Random'=clo(c(tbR)),'PVA'=clo(c(tbMV)),'F-stat'=clo(c(tbF)),'GPF'=clo(c(tbGPF)),
                  'GraphBM'=clo(c(tbBM)),'Fisher'=clo(c(tbFish)),'Voronoi'=clo(c(tbVoronoi)))

melt(ddf,"Voronoi") %>%
ggplot(aes(Voronoi, value, colour = variable)) + 
  geom_point(size=3) +
  facet_wrap(~ variable, nrow=1)+
  geom_abline(slope = 1,intercept = 0)+
  scale_y_continuous(limits = c(0.002,0.03))
ggsave('Voronoi_vs_Obj_Fcn.tiff',height = 3,width=16)

# save(list=ls(),file='graph_partitioning_analysis_big_workspace_final')














# Bias among unevenly sized subtrees --------------------------------------

### We'll partition tree into a set of 5 partitions ranging from
### 19:21, 12:28, 6:35, and 3:37
### corresponding to edges 23, 24, 10 and 3, respectively
### and nodes 53, 54, 47, and 56, respectively
Vmax2 <- function(reps,V){
  A <- matrix(rnorm(reps*ncol(V)),nrow=ncol(V))
  mx <- (V %*% A) %>% apply(.,MARGIN=2,FUN=function(x) which.max(abs(x)))
  tb <- table(mx)
  return(tb)
}

gpf2 <- function(reps,Grps,tree,samples=as.character(1:10),prob=0.5){
  mx <- numeric(reps)
  n <- length(samples)
  spp <- tree$tip.label
  D <- length(spp)
  SampleFrame  <- expand.grid(samples,c('R','S'))
  names(SampleFrame) <- c('sample','G')
  SampleFrame <- as.data.table(SampleFrame)
  for (rr in 1:reps){
    PA <- lapply(spp,FUN=function(x,samples,prob) samples[as.logical(rbinom(n,1,prob))],samples,prob)
    DF <- data.table('Species'=rep(spp,sapply(PA,length)),'sample'=unlist(PA),'N'=1,stringsAsFactors = F)
    X <- data.table('sample'=samples,'zz'=rnorm(n))
    mx[rr] <- which.max(sapply(Grps,getObjective,tree,DF,X,SampleFrame,frmla=cbind(Successes,Failures)~zz*G))
  }
  
  tb <- table(mx)
  return(tb)
}

fishtest2 <- function(reps,Grps,tree,p=0.5){
  Ftest <- function(grps,Z){
    s1 <- Z[grps[[1]]]
    s2 <- Z[grps[[2]]]
    n1 <- sum(s1)
    n2 <- sum(s2)
    ft <- fisher.test(matrix(c(n1,length(s1)-n1,n2,length(s2)-n2),ncol=2),alternative='two.sided')
    return(log(ft$p.value))
  }
  TestFunction <- function(Z,Grps){
    return(sapply(Grps,Ftest,Z=Z))
  }
  n <- Ntip(tree)
  E <- numeric(reps)
  
  for (i in 1:reps){
    Z <- rbinom(n,1,p)
    E[i] <- TestFunction(Z,Grps) %>% which.min
  }
  return(table(E))
}

tableListCount <- function(tb,ix1,ix2){
  tbcnt <- function(tb,ix1,ix2){
    n1 <- sum(tb[as.character(ix1)],na.rm = T)
    return(n1)
  }
  cnt <- sapply(tb,tbcnt,ix1,ix2) %>% sum
  return(cnt)
}

### we'll do three phylofactorizations with null data: (1) argmax(v'X), (2) gpf (3) FisherFactor

### We want to know the probability of drawing one tree or another relative to the number of edges in each tree
nodes <- c(53,54,47,56)
methods <- c('Voronoi','gpf','Fisher')
reps <- rep(3000,7)
spp <- tree$tip.label
cl <- phyloFcluster(7)

TwoTreeData <- data.table('node'=rep(nodes,each=3),'method'=rep(methods,4),
                          'N1'=numeric(12),'N2'=numeric(12),'Nedges1'=numeric(12),'Nedges2'=numeric(12),stringsAsFactors = F)

for (nn in nodes){
  g1 <- phangorn::Descendants(tree,nn,'tips')[[1]]
  
  t1 <- drop.tip(tree,setdiff(tree$tip.label,tree$tip.label[g1]))
  t2 <- drop.tip(tree,tree$tip.label[g1])
  Grps1 <- getPhyloGroups(t1)
  Grps1 <- Grps1[2:length(Grps1)]
  Grps2 <- getPhyloGroups(t2)
  Grps2 <- Grps2[2:length(Grps2)]
  ix1 <- 1:length(Grps1)
  ix2 <- 1:length(Grps2)+length(Grps1)
  Grps <- c(Grps1,Grps2)
  
  TwoTreeData$Nedges1[TwoTreeData$node==nn] <- length(ix1)
  TwoTreeData$Nedges2[TwoTreeData$node==nn] <- length(ix2)
  
  V <- sapply(Grps,ilrvec,D) %>% t
  TBV <- parLapply(cl,reps,Vmax2,V)
  TBgpf <- parLapply(cl,reps,gpf2,Grps,tree)
  TBF <- parLapply(cl,reps,fishtest2,Grps,tree)
  
  TwoTreeData[node==nn & method=='Voronoi','N1'] <- tableListCount(TBV,ix1,ix2)
  TwoTreeData[node==nn & method=='gpf','N1'] <- tableListCount(TBgpf,ix1,ix2)
  TwoTreeData[node==nn & method=='Fisher','N1'] <- tableListCount(TBF,ix1,ix2)
  
}

TwoTreeData$N2 <- sum(unlist(reps))-TwoTreeData$N1
TwoTreeData[,RelativeSize:=Nedges1/Nedges2]
TwoTreeData[,Score:=(N1/N2)/RelativeSize]
TwoTreeData[,Unevenness:=1-RelativeSize]
TwoTreeData[,method:=factor(method)]

stopCluster(cl)
rm('cl')


ggplot(TwoTreeData,aes(x=Unevenness,y=Score,group=method,color=method))+
  geom_line(lwd=2)+
  geom_point(cex=4)+
  geom_abline(intercept=1,slope=0,lwd=2)+
  scale_y_continuous(name='Small subtree bias')+
  scale_x_continuous(name='1-r/s')
  ggtitle('Relative bias for drawing small sub-graph')
ggsave('Bias_btwn_subgraphs.tiff',height=3,width=3)


save(list=ls(),file='graph_partitioning_workspace_FINAL2')
