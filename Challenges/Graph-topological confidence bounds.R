######### graph confidence intervals
library(phylofactor)
library(phangorn)
library(parallel)
setwd('Challenges')

set.seed(1)
D=30
n=10  ## will use two different values of n
tree <- rtree(D)

nd=42
grp1 <- Descendants(tree,nd,'tips')[[1]]
grp2 <- setdiff(1:D,grp1)

Grps <- getPhyloGroups(tree)
Grps <- Grps[2:length(Grps)]  #omit redundant first group
V <- sapply(Grps,ilrvec,D)
edgs <- apply(V,2,getFactoredEdges,tree=tree)
root.edges <- edgs[sapply(edgs,length)>1][[1]]
edgs[sapply(edgs,length)>1]=root.edges[2]  #our reference edge for the root.
V <- t(V)

pf.rep <- function(reps,tree,grp1,grp2,n,V,edgs,dz=0.05){
  
  winners <- numeric(reps)
  for (rr in 1:reps){
    x <- rnorm(n)
    mu1 <- (dz/2)*x
    mu2 <- -(dz/2)*x
    Y <- matrix(NA,nrow=Ntip(tree),ncol=n)
    for (i in 1:nrow(Y)){
      if (i %in% grp1){
        Y[i,] <- rnorm(n,mean=mu1)
      } else {
        Y[i,] <-  rnorm(n,mean=mu2)
      }
    } 

    Y <- V %*% Y
    GG <- apply(Y,1,FUN=function(y,x) glm(y~x),x=x)
    winners[rr] <- edgs[which.max(sapply(GG,FUN=function(gg) getStats(gg)['F']))]
  }
  return(winners)
}

reps <- rep(1e5,7)
cl <- phyloFcluster(7)
Winners <- parSapply(cl,reps,pf.rep,tree=tree,grp1=grp1,grp2=grp2,n=n,V=V,edgs=edgs)
Winners1 <- parSapply(cl,reps,pf.rep,tree=tree,grp1=grp1,grp2=grp2,n=n,V=V,edgs=edgs,dz=0.1)
Winners2 <- parSapply(cl,reps,pf.rep,tree=tree,grp1=grp1,grp2=grp2,n=n,V=V,edgs=edgs,dz=0.2)
Winners5 <- parSapply(cl,reps,pf.rep,tree=tree,grp1=grp1,grp2=grp2,n=n,V=V,edgs=edgs,dz=0.5)
stopCluster(cl)
rm('cl')

save(list=ls(),file='graph_topology_of_errors')


# Visualization of graph-topology of error --------------------------------
# load('graph_topology_of_errors')
true.edg <- which(tree$edge[,2]==nd)

### counting & sorting edges
finishtbl <- function(Winners,edgs){

  tb <- table(unlist(Winners))
  
  n0 <- sum(!(sapply(edgs,toString) %in% names(tb)))
  tb0 <- rep(0,n0)
  names(tb0) <- unlist(edgs)[!(sapply(edgs,toString) %in% names(tb))]
  
  tb <- c(tb,tb0)
  tb <- tb[order(as.numeric(names(tb)))]
  return(tb)
} 

getPscale <- function(tb){

  probs <- tb/sum(tb)
  probs <- c('1'=probs['16'],probs)
  names(probs)[1] <- '1'
  
  pscale <- probs/max(probs)
  return(pscale)
}
getProb <- function(tb){
  
  probs <- tb/sum(tb)
  probs <- c('1'=probs['16'],probs)
  names(probs)[1] <- '1'
  return(probs)
}

probs <- finishtbl(Winners,edgs) %>% getProb
probs1 <- finishtbl(Winners1,edgs) %>% getProb
probs2 <- finishtbl(Winners2,edgs) %>% getProb
probs5 <- finishtbl(Winners5,edgs) %>% getProb

pscale <- finishtbl(Winners,edgs) %>% getPscale
pscale1 <- finishtbl(Winners1,edgs) %>% getPscale
pscale2 <- finishtbl(Winners2,edgs) %>% getPscale
pscale5 <- finishtbl(Winners5,edgs) %>% getPscale

cols <- rgb(red = pscale,green=0,blue = 1-pscale)
cols1 <- rgb(red = pscale1,green=0,blue = 1-pscale1)
cols2 <- rgb(red = pscale2,green=0,blue = 1-pscale1)
cols5 <- rgb(red = pscale5,green=0,blue = 1-pscale5)

tiff('Graph_topology_error.tiff',height=600,width=1600)
  par(mfrow=c(2,3))
  # plot.phylo(tree,edge.color = cols,edge.width=2+4*pscale)
  plot.phylo(tree,edge.color = cols1,edge.width=2+4*pscale1,main='effect size 0.1')
  edgelabels('*',true.edg,cex=4)
  plot.phylo(tree,edge.color = cols2,edge.width=2+4*pscale2,main='effect size 0.2')
  edgelabels('*',true.edg,cex=4)
  plot.phylo(tree,edge.color = cols5,edge.width=2+4*pscale5,main='effect size 0.5')
  edgelabels('*',true.edg,cex=4)
  
  
  
  ## true.edg-1 ID's the correct edge - plot distances vs. pscale
  dsts <- acos(abs((V %*% t(V))[true.edg-1,]))
  ix <- 2:length(pscale)
  plot(dsts,probs1[ix],pch=16,cex=3,cex.axis=2,log='y',main='effect size 0.1',xlab='Angular distance from e*',ylab='P{edge}')
  plot(dsts,probs2[ix],pch=16,cex=3,cex.axis=2,log='y',main='effect size 0.1',xlab='Angular distance from e*',ylab='')
  plot(dsts,probs5[ix],pch=16,cex=3,cex.axis=2,log='y',main='effect size 0.1',xlab='Angular distance from e*',ylab='')

dev.off()