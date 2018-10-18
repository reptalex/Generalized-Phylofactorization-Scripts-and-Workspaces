library(parallel)
library(phylofactor)
library(ape)
library(phangorn)

setwd('Challenges/')
source('edge_partitioning_functions.R')

####### testing voronoi-cell-approximation for probabilities of P(edge) vs. unevenness

## three key variables:
# 1) D: tree size
# 2) nn: number of nearest-neighbors over which we average angles to approximate voronoi cell size
#          - optional weighted average, e.g. for kth NN, w(k)=1/k. 
# 3) a: exponent we raise the angles to minimize chisq.test(table(edges),p=angles^a)

## to start off, we will use the maximum-variance objective function, i.e. PhyCA
## Then, we'll consider fastBM

Dset <- c(10,20,30) #,40,80) #,160,320)
nnset <- as.list(1:16)
aset <- seq(0,10,by=0.01) ## will find ahat by two ways:
 #1) b=mean(tbMV[sapply(tips,toString)])/mean(tbMV[sapply(setdiff(2:(E+1),tips),toString)])
 #   ahat=a[which.min(abs(b-alphamn(a,tips,E)))]
 #2) ahat=a[which.max(chisq.test(edgs,p=clo(angles^a)))]

#1)
min.t2b <- function(ms,aset,tips,E,t2b,tbl){
  clo <- function(x) x/sum(x)
  obs <- intersect(names(tbl),2:E) %>% sapply(.,toString)
  ts <- intersect(tips,obs)
  bs <- setdiff(obs,ts)
  
  a=aset[which.min(abs(t2b-alphamn(aset,tips,E,clo(ms))))]
  P <- chisq.test(tbl[obs],p=clo(ms[obs]^a))$p.value
  return(list('max'=a,'chisqPval'=P))
}

#2)
min.chisq <- function(ms,aset,tbl,E,tips){
  clo <- function(x) x/sum(x)
  obs <- intersect(names(tbl),2:E) %>% sapply(.,toString)
  ts <- intersect(tips,obs)
  bs <- setdiff(obs,ts)
  ch <- sapply(aset,FUN=function(a,ms,tbl) chisq.test(tbl,p=clo(ms^a))$p.value,ms=ms[obs],tbl=tbl[obs])
  return(list('max'=aset[which.max(ch)],'chisqPval'=max(ch)))  
}

n=2000                         #number of replicate phylofactorizations per core
ncores=7                       #number of cores
reps <- as.list(rep(n,ncores))
cl <- phyloFcluster(ncores)
clusterEvalQ(cl,expr = source('edge_partitioning_functions.R'))

EdgTrees <- vector(mode='list',length=length(Dset))
names(EdgTrees) <- sapply(Dset,FUN=function(D) paste('D=',toString(D),sep=''))

# PhyCA -------------------------------------------------------------------


for (D in Dset){
  d <- match(D,Dset)
  set.seed(1)
  tree=rtree(D)
  E <- Nedge(tree)
  tree$edge.length=rep(0.1,E)
  tips <- which(tree$edge[,2] %in% 1:D) %>% sapply(.,toString)
  edgs <- parLapply(cl,X=reps,fun=pf.edges,tree=tree,choice='var',method='max.var')
  edgs <- unlist(edgs)
  edgMV <- edgs[!edgs==1]
  
  EdgTrees[[d]]$EdgeTable <- table(edgMV)
  EdgTrees[[d]]$tree <- tree
  EdgTrees[[d]]$t2bProb <- tip_to_basal_prob(edgMV,tree)
  
  angles <- getAngles(tree)
  
  nms <- sapply(nnset,FUN=function(nn) paste('nn=',toString(nn),sep=''))
  means <- nnset %>% sapply(.,FUN=function(nn,A) angleMean(A,nn,weighted=F),A=angles)
  Wmeans <- nnset %>% sapply(.,FUN=function(nn,A) angleMean(A,nn,weighted=T),A=angles)
  colnames(means) <- nms
  colnames(Wmeans) <- nms
  
  EdgTrees[[d]]$angleMeans <- means
  EdgTrees[[d]]$WangleMeans <- Wmeans
  
  ### for each nn, we find the min ahat using either method 1 or 2.
  t2b <- tip_to_basal_prob(edgs,tree) ## tip-to-basal prob.
  
  for (w in c(F,T)){
    if (w){
      T2B <- suppressWarnings(apply(Wmeans,2,min.t2b,aset=aset,tips=tips,E=E,t2b=t2b,tbl=EdgTrees[[d]]$EdgeTable))
      CHSQ <- suppressWarnings(apply(Wmeans,2,min.chisq,aset=aset,tbl=EdgTrees[[d]]$EdgeTable,E=E,tips=tips))
    } else {
      T2B <- suppressWarnings(apply(means,2,min.t2b,aset=aset,tips=tips,E=E,t2b=t2b,tbl=EdgTrees[[d]]$EdgeTable))
      CHSQ <- suppressWarnings(apply(means,2,min.chisq,aset=aset,tbl=EdgTrees[[d]]$EdgeTable,E=E,tips=tips))
    }
    dd <- data.frame('nn'=unlist(nnset),
                     'a'=sapply(T2B,FUN=function(x) x$max),
                     'P'=sapply(T2B,FUN=function(x) x$chisqPval),
                     'method'=rep('T2B',length(nnset)))
    dd2 <- data.frame('nn'=unlist(nnset),
                      'a'=sapply(CHSQ,FUN=function(x) x$max),
                      'P'=sapply(CHSQ,FUN=function(x) x$chisqPval),
                      'method'=rep('CHSQ',length(nnset)))
    if (w){
      EdgTrees[[d]]$weighted_nn_a <- rbind(dd,dd2)
    } else {
      EdgTrees[[d]]$unweighted_nn_a <- rbind(dd,dd2)
    }
  }
  
}

stopCluster(cl)
rm('cl')

save(list=ls(),file='Voronoi_Workspace')



# Visualizing results -----------------------------------------------------

load(file='Voronoi_Workspace')

names(EdgTrees[[1]])
names(EdgTrees[[1]]$unweighted_nn_a)
### first, what the maximum t2b

sapply(EdgTrees,FUN=function(d) d$t2b)
## the t2b seems to decrease, but this may be due to the sparse sampling

## how did weighted vs. unweighted do?
sapply(EdgTrees,FUN=function(d) d$unweighted_nn_a[which.max(d$unweighted_nn_a$P),]) %>% t
sapply(EdgTrees,FUN=function(d) d$weighted_nn_a[which.max(d$weighted_nn_a$P),]) %>% t

### no clear pattern except increasing n, decreasing a and high P-values
######### need a very very large sample size ot get a handle on this.

#### probably want something like 90,000- want to get good coverage of all edges.
#### To do this formally, we need to think: what is the P{seen_edges==E|N}
#### and get P{seen_edges==E|n} to be near 1. 

#### Let's compute this by induction on P{n+1 is unique|n,obs}
#### E[time to next unique] = (no. obs edges)/(no. edges)
#### e.g. E[t|n=1]=1
####      E[t|obs=n]=E/(E-n)  - if 300 edges 30 obs, p=0.9 of unique and so 1/0.9=10/9
####
fn <- function(x,E) E/(E-x)
sum(fn(1:(E-1),E))   ###### 4488. Let's sample on order of 14,000 or 2,000 reps per core.




# Vcell estimation -------------------------------------------------
### let's see if edges are drawn at a rate equal to the
### probability of the edge having the maximum projection of a random vector
load('Voronoi_Workspace')

rproj <- function(n,V,output.max=T){
  M <- matrix(rnorm(nrow(V)*n),ncol=nrow(V))
  proj <- abs(M %*% V)
  if (output.max){
    maxs <- apply(proj,1,which.max)
    return(maxs)
  } else {
    return(proj)
  }
}

n=1e5
ncores=7
reps=as.list(rep(n,ncores))
cl <- makeCluster(ncores)
clusterExport(cl,varlist='rproj')

par(mfrow=c(1,3))
for (D in Dset){
  
  ix <- match(D,Dset)
  tree <- EdgTrees[[ix]]$tree
  grps <- getPhyloGroups(tree)
  V <- sapply(grps,ilrvec,n=D) #rows are species, columns are groups
  V <- V[,2:(2*D-2)]  # drops redunant root
  
  colnames(V) <- 2:(2*D-2)
  rownames(V) <- tree$tip.label
  
  
  
  ########## parallelizing computation
  
  maxs <- parLapply(cl,reps,rproj,V=V) %>% unlist
  tb <- table(maxs)
  names(tb) <- 2:(2*D-2)
  tb <- tb[names(EdgTrees[[ix]]$EdgeTable)]
  
  ### double check that we get 1-1 for replicates
  # tb2 <- parLapply(cl,reps,rproj,V=V) %>% unlist %>% table
  # plot(c(tb),c(tb2))
  nn <- length(tb)
  plot(c(tb)/sum(tb),c(EdgTrees[[ix]]$EdgeTable)/sum(EdgTrees[[ix]]$EdgeTable),main=paste('D=',toString(D),sep=''))
  abline(0,1,lwd=2,col='blue')
}

stopCluster(cl)
rm('cl')

save(list=ls(),file='Voronoi_Workspace')