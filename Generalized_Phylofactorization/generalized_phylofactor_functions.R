
logit <- function(p) log(p/(1-p))
ilogit <- function(x) exp(x)/(1+exp(x))
makeContrastFrame <- function(grp,Y,Z,tree,igHomogeneity=T){
  R <- tree$tip.label[grp[[1]]]
  S <- tree$tip.label[grp[[2]]]
  
  if (igHomogeneity){
    r=length(R)
    s=length(S)
    YR <- colSums(Y[R,,drop=F])
    YS <- colSums(Y[S,,drop=F])
    DF <- data.frame('successes'=c(YR,YS),
                     'failures'=rep(c(r-YR,s-YS),times=length(YR)),
                     'group'=rep(c('R','S'),each=length(YR)),
                     'sample'=rep(colnames(Y),2))
  } else {
    DF <- data.frame('sample'=rep(colnames(Y),times=nrow(Y)),
                     'successes'=c(t(Y)),
                     'failures'=1-c(Y),
                     'id'=rep(rownames(Y),each=ncol(Y)))
    DF$group <- rep('S',nrow(DF))
    DF$group[DF$id %in% R] <- 'R'
    DF$group <- factor(DF$group)
  }
  
  DF <- cbind(DF,Z[match(DF$sample,Z$sample),])
  return(DF)
}
binomialContrast <- function(grp,Y,Z,tree,igHomogeneity=T){
  DF <- makeContrastFrame(grp,Y,Z,tree,igHomogeneity)
  ContrastGlm <- glm(cbind(successes,failures)~x+group*y+group*z,family=binomial(link='logit'),data=DF)
  return(ContrastGlm)
}

FactorContrast <- function(grp,Y,Z,tree,igHomogeneity=T){
  cGlm <- binomialContrast(grp,Y,Z,tree,igHomogeneity)
  omega <- anova(cGlm,test='Chisq')['group:z','Deviance']
  return(omega)
}


randomFactorSim <- function(reps,dz,sd.ig,tree,Grps,n){
  clid <- signif(rnorm(1),digits = 3)
  OUTPUT <- data.frame('true_group'=numeric(2*reps),
                       'factored_group'=numeric(2*reps),
                       'correct_group'=logical(2*reps),
                       'igHomogeneity'=rep(c(T,F),each=reps),
                       'dz_true'=rep(dz,2*reps),
                       'dz_est'=numeric(2*reps),
                       'sd'=rep(sd.ig,2*reps),
                       'sim'=rep(1:reps,times=2),
                       'cluster_id'=rep(clid,2*reps))
  
  nGrps <- length(Grps)
  D <- Ntip(tree)
  Y <- matrix(NA,nrow=D,ncol=n)
  rownames(Y) <- tree$tip.label
  colnames(Y) <- 1:n
  for (rr in 1:reps){
    ix <- sample(nGrps,1)
    OUTPUT$true_group[c(rr,reps+rr)]=ix
    grp1 <- tree$tip.label[Grps[[ix]][[1]]]
    grp2 <- tree$tip.label[Grps[[ix]][[2]]]
    
    Z <- data.frame('sample'=1:n,'x'=rnorm(n,sd=0.3),'y'=rnorm(n,sd=0.3),'z'=rnorm(n))
    etaR <- Z$x+Z$y+(0.1+dz/2)*Z$z
    etaS <- Z$x-Z$y+(0.1-dz/2)*Z$z
    
    
    for (i in tree$tip.label){
      if (i %in% grp1){
        Y[i,] <- rbinom(n=n,size=1,prob=ilogit(etaR+rnorm(n,sd=sd.ig)))
      } else {
        Y[i,] <- rbinom(n=n,size=1,prob=ilogit(etaS+rnorm(n,sd=sd.ig)))
      }
    }
    weirdos <- (rowSums(Y) %in% c(0,n))
    while (any(weirdos)){
      for (ii in rownames(Y)[weirdos]){
        if (ii %in% grp1){
          Y[ii,] <- rbinom(n=n,size=1,prob=ilogit(etaR+rnorm(n,sd=sd.ig)))
        } else {
          Y[ii,] <- rbinom(n=n,size=1,prob=ilogit(etaS+rnorm(n,sd=sd.ig)))
        }
      }
      weirdos <- (rowSums(Y) %in% c(0,n))
    }
    
    
    Omegas <- sapply(Grps,FUN=FactorContrast,Y=Y,Z=Z,tree=tree,igHomogeneity=T)
    winner=which.max(Omegas)
    OUTPUT$factored_group[rr]=winner
    gg=binomialContrast(Grps[[winner]],Y,Z,tree)
    OUTPUT$dz_est[rr] <- coef(gg)['groupS:z']
    
    
    Omegas <- sapply(Grps,FUN=FactorContrast,Y=Y,Z=Z,tree=tree,igHomogeneity=F)
    winner=which.max(Omegas)
    OUTPUT$factored_group[rr+reps]=winner
    gg=binomialContrast(Grps[[winner]],Y,Z,tree,igHomogeneity=F)
    OUTPUT$dz_est[rr+reps] <- coef(gg)['groupS:z']
    
  }
  OUTPUT$correct_group <- OUTPUT$true_group==OUTPUT$factored_group
  
  return(OUTPUT)
}