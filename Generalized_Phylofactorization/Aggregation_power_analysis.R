# Power Analysis ----------------------------------------------------------

## in this section, we simulate many affect sizes, with varying degrees of intragroup homogeneity.
## the main function produces an effect along a random edge for many replicates.
## the output will be a data frame: 

## OUTPUT = true_group, factored_group, correct_group, igHomogeneity, dz_true, dz_est, sim

library(phylofactor)
library(parallel)
setwd('Generalized_Phylofactorization')
source('generalized_phylofactor_functions.R') 
##these functions in the sourced file & in the script below became the scaffold
## for the current gpf() function
set.seed(1)
D=15
tree <- rtree(D)
n=60
DZ <- seq(0,3,length.out=9) ### effect sizes
SD <- c(0,1,2)              ### intra-group homogeneity

Grps <- getPhyloGroups(tree)     ### obtain groups just once - each item contains groups R-S for each edge.
mrnode <- function(grp,tree){   #gets ancestral node for each edge
  grp <- grp[[1]]
  if (length(grp)==1){
    output=grp
  } else {
    output <- getMRCA(tree,tree$tip.label[grp])
  }
  return(output)
}
Grp_nodes <- sapply(Grps,mrnode,tree=tree)

cl <- makeCluster(7)
clusterExport(cl,varlist=c('ilogit','makeContrastFrame','binomialContrast','FactorContrast','randomFactorSim'))
clusterEvalQ(cl,library(ape))


reps <- rep(100,7)
DF <- NULL
for (dz in DZ){
  for (sd.ig in SD){
    
    
    OUTPUT <- parLapply(cl,reps,fun=randomFactorSim,dz=dz,sd.ig=sd.ig,tree=tree,Grps=Grps,n=n)
    
    if (is.null(DF)){
      DF <- OUTPUT[[1]]
      for (i in 2:length(OUTPUT)){
        DF <- rbind(DF,OUTPUT[[i]])
      }
    } else {
      for (i in 1:length(OUTPUT)){
        DF <- rbind(DF,OUTPUT[[i]])
      }
    }
    
  }
}

np <- function(ix,Grp_nodes,tree){
  n1 <- Grp_nodes[ix[1]]
  n2 <- Grp_nodes[ix[2]]
  if (n1==n2){
    return(NULL)
  } else {
    output <- nodepath(tree,n1,n2)
    return(output)
  }
}

nodepaths <- apply(DF[,c('true_group','factored_group')],MARGIN=1,FUN=np,Grp_nodes=Grp_nodes,tree=tree)
nodepath_lengths <- sapply(nodepaths,length)
plot(ecdf(nodepath_lengths))

DF$EdgeDist <- nodepath_lengths

rnodes <- sample(length(Grps),size = nrow(DF)*2,replace=T) %>% matrix(.,ncol=2)
rnps <- apply(rnodes,MARGIN=1,np,Grp_nodes=Grp_nodes,tree=tree)
rnpsL <- sapply(rnps,length)
plot(ecdf(rnpsL))
lines(ecdf(nodepath_lengths),col='blue')


### forgot to add "sd.ig" to the data frame - this will do it
## for each DZ, we cycle through SD. 
## Then, for each (DZ,SD), we cycle through 700 replicates for each of two igHomogeneity
## totalling 1400 rows for each (DZ,SD). These are captured in rep(rep(SD,each=2),each=700) = ( (each igHomogeneity) each replicate)
## Then, we repeat that pattern a number of times  (... times length(DZ))
# sd.df <- rep(rep(rep(SD,each=2),each=sum(reps)),times=length(DZ))

# DF$sd <- sd.df
DF$n <- rep(n,nrow(DF))

stopCluster(cl)
rm('cl')
save(list=ls(),file=paste('Bernoulli_aggregation_workspace_n',toString(n),sep=''))



############# summarizing results:

X <- NULL
for (n in c(5, 10, 30, 60)){
  load(paste('Bernoulli_aggregation_workspace_n',toString(n),sep=''),envir = )
  if (n==5){
    X <- DF
  } else {
    X <- rbind(X,DF)
  }
}

DF <- X
rm('X')

DF$n <- as.factor(DF$n)
levels(DF$n) <- c('n=5','n=10','n=30','n=60')

DF$igHomogeneity <- as.factor(DF$igHomogeneity)
levels(DF$igHomogeneity) <- c('Bernoulli','Binomial Aggregation')

require(reshape2)
require(ggplot2)
df_melt <- melt(DF, id = c("dz_true","igHomogeneity","sd","n"))
df <- dcast(df_melt, igHomogeneity + dz_true + sd +n ~ variable, mean)
df$cg_numeric = df$correct_group

ggplot(df,aes(x=dz_true,y=correct_group,group=sd,color=sd))+
  geom_line()+
  geom_smooth()+
  facet_wrap(~n+igHomogeneity,nrow=4)

DF$cg_numeric <- as.numeric(DF$correct_group)

gg1 <- ggplot(DF,aes(x=dz_true,y=cg_numeric,group=sd,color=sd))+
  geom_line(data=df,lwd=2)+
  # geom_smooth(data=DF,method='gam',method.args=list(family=binomial))+
  facet_wrap(~n+igHomogeneity,nrow=4)+
  scale_x_continuous(name='Effect Size')+
  scale_y_continuous(name='P{Correct edge}')
gg1
ggsave('Binomial_aggregation_Pedge.png',width=12,height=8)


gg2 <- ggplot(DF,aes(x=dz_true,y=EdgeDist,group=sd,color=sd))+
  geom_line(data=df,lwd=2)+
  # geom_smooth(data=DF,method='gam',method.args=list(family=poisson))+
  facet_wrap(~n+igHomogeneity,nrow=4)+
  scale_x_continuous(name='Effect Size')+
  scale_y_continuous(name='Distance to Correct Edge')
gg2
ggsave('Binomial_aggregation_DistToEdge.png',width=12,height=8)

tiff('Binomial_Bernoulli_Performance_all.tiff',height=1200,width=1000)
gridExtra::grid.arrange(gg1,gg2,layout_matrix=matrix(c(1,2),ncol=2))
dev.off()
# ggsave('Binomial_aggregation_performance.png',width=12,height=8)

### Estimates of dz:
ggplot(DF,aes(x=dz_true,y=dz_est,group=sd,color=sd))+
  geom_boxplot()+
  facet_grid(~n+igHomogeneity,nrow=4)
