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

edge.chisqtest <- function(edg,...){
  tb <- table(edg)
  ct <- chisq.test(tb,...)
  return(ct$p.value)
}

tip_to_basal_prob <- function(edg,tree){
  tips <- which(tree$edge[,2] %in% 1:Ntip(tree)) %>% sapply(.,toString)
  tbl <- table(edg[!edg==1])
  
  #restrict ourselves to observed edges
  obs <- intersect(2:(Nedge(tree)),names(tbl))
  ts <- intersect(tips,obs)
  bs <- setdiff(obs,ts)
  tip.mean=mean(tbl[ts],na.rm=T)
  basal.mean=mean(tbl[bs],na.rm=T)
  return(tip.mean/basal.mean)
}

alphamn <- function(a,tips,E,ms){
  sapply(a,FUN=function(a,tips,E) mean(ms[tips]^a)/mean(ms[setdiff(1:E,tips)]^a),tips=tips,E=E)
}


getAngles <- function(tree){
  V <- phylofactor::getPhyloGroups(tree) %>% sapply(.,FUN=function(x,D) phylofactor::ilrvec(x,D),D=D)
  colnames(V) <- 1:ncol(V)
  V <- V[,2:ncol(V)]   ### omits the one redundant root tip.
  projV <- t(V) %*% V
  projV[abs(projV)>1]=1
  angles <- acos(projV)   ### 
  diag(angles) <- 0
  
  angles <- -abs(angles-pi/2)   ### this is a measure of similarity
  angles <- angles-min(angles)
  angles <- angles*(180/pi)
  colnames(angles) <- colnames(V)
  rownames(angles) <- colnames(V)
  return(angles)
}

angleMean <- function(angles,nn,weighted=F){
  return(apply(angles,1,FUN=function(x) weighted.mean(sort(x[x>0],decreasing = F)[1:nn],w = (1/1:nn)^weighted)))
}
