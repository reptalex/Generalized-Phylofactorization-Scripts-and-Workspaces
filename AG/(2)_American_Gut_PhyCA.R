######### American Gut Phylofactorization
library(phylofactor)
library(xlsx)
library(ggplot2)
library(mgcv)
library(ggtree)

setwd('AG')
load('AmericanGut_Workspace')

### finally we'll focus on feces
ix <- which(X$BODY_SITE=='UBERON:feces')
X <- X[ix,]
OTUTable <- OTUTable[,sapply(X$SampleID,toString)]
rm('ix')

### We'll focus on more accurate samples & OTUs
### We'll trim entries to those with >50,000 sequence counts
OTUTable <- OTUTable[,colSums(OTUTable)>5e4]
OTUTable <- OTUTable[rowSums(OTUTable)>=ncol(OTUTable),]
tree <- drop.tip(tree,setdiff(tree$tip.label,rownames(OTUTable)))
OTUTable <- OTUTable[tree$tip.label,]


phca <- PhyCA(OTUTable,tree,ncores=7,ncomponents = 10)

save(list=ls(),file='American_gut_PhyCA')



# Visualization and analysis ----------------------------------------------
phca$factors

### the following uses the minimu set of shortest unique taxonomic prefixes to summarize groups
n=1
s <- summary(phca,Taxonomy,n)  ## pf.tidy returns a minimal list of unique taxonomic identifiers
s$taxa.split$group1  ### taxonomic detail of group

### Doing so yields the following:
# n=1    #1229-member monophyletic clade of Firmicutes c_Erysipelotrichi etc.
# n=2    #215-member monophyletic clade of f__Lachnospiraceae partitioned from previous clade
# n=3    #81 Bacteroidetes of genus Bacteroides
# n=4    #60 Ruminococcaceae
# n=5    #27 members of f_[tissierellaceae]
# n=6    #5 member clade of g_Prevotella s_copri
# n=7    #41 Gammaproteobacteria o_Enterobacteriales
# n=8    #15 Actinobacteria f_Bifidobacteriaceae (g_Bifidobacterium & g_Gardnerella)
# n=9    #32 o_Clostridiales, unclassified family
# n=10   #64 o_Clostridiales, unclassified genera

## which we simplify below:
FactorNames <- c('Firmicutes_Erysipelotrichi_etc',
                 'Lachnospiraceae_genera',
                 'Bacteroides_species',
                 'Ruminococcaceae',
                 'Tissierellaceae',
                 'Prevotella_copri',
                 'Gammaproteobacteria_Enterobacteriales',
                 'Bifidobacteriaceae_genera',
                 'Clostridiales_f_',
                 'Clostridiales_g_')



#### this needs to be done for ggtree plotting.
phca$tree$edge.length[is.na(phca$tree$edge.length)] <- min(phca$tree$edge.length[!is.na(phca$tree$edge.length)])
## errors in the tree edge.lengths require this step for a visualization of the tree

#### phylogeny plot done with pf.tree
ggtr <- pf.tree(phca,bg.color = NA)
lgnd <- ggtr$legend
ggtr$ggplot
ggsave('American_Gut_PhyCA_transparent.png',height=8,width=8)
tiff('American_Gut_PhyCA_legend.tiff',height=800,width=800)
  plot(rep(1,10),10:1,col=lgnd$colors,pch=16,cex=8,xlim=c(1,5))
  legend(2,10,legend = c('Firms','Lachnos','Bacteroides','Ruminococs','tissierell','prevotella','Enterobacteriales','Bifido','f_','g_'),
         col = lgnd$colors,pch=16,cex=2.2)
dev.off()

##### Limited phylogeny plot of clades hilighted in paper
ggtr <- pf.tree(phca,factors = c(1,6,7),bg.color = 'grey',rotation=90)
cls <- ggtr$legend$colors
rotate_tree(ggtr$ggplot,-90)
ggsave('American_Gut_PhyCA2.png',height=8,width=8)
tiff('American_Gut_PhyCA_legend2.tiff',height=800,width=800)
plot(rep(1,3),3:1,col=lgnd$colors,pch=16,cex=8,xlim=c(1,5))
legend(2,3,legend = c('Firms','prevotella','Enterobacteriales'),
       col = cls,pch=16,cex=2.2)
dev.off()
########################################################







# Regression on component scores ------------------------------------------
### align sample IDs of X with columns of OTUtable
X <- X[match(colnames(OTUTable),X$SampleID),]


############ Creating explanatory variables #############
plants <- X$TYPES_OF_PLANTS %>% relevel(.,'6 to 10') %>% relevel(.,'Less than 5')
age <- as.numeric(sapply(X$AGE_YEARS,toString))
bmi <- as.numeric(X$BMI)
alcohol <- X$ALCOHOL_FREQUENCY %>% 
  relevel(.,'Never') %>%
  relevel(.,'Rarely (a few times/month)') %>%
  relevel(.,'Regularly (3-5 times/week)') %>%
  relevel(.,'Occasionally (1-2 times/week)') %>%
  relevel(.,'Daily')
sex <- X$SEX
sex[!sex%in%c('female','male')] <- NA
sex <- factor(sex)

ABX <- sapply(X$ANTIBIOTIC_HISTORY,toString)
ABX[ABX=='I have not taken antibiotics in the past year.'] <- 'Longer'
ABX <- factor(ABX,levels = c('Unspecified','Longer','Year','6 months','Month','Week'))
ibd <- X$SUBSET_IBD
ibd[ibd=='False'] <- 'false'
ibd[ibd=='True'] <- 'true'
ibd <- factor(ibd)
###########################

########################### useful functions #######################
factorILR <- function(nn){
  y <- t(phca$basis[,nn]) %*% log(phca$Data) %>% t
  return(y)
}
getReg <- function(nn,model=F){
  y <- factorILR(nn)
  gg <- glm(y~plants+age+bmi+alcohol+sex+ABX+ibd)  ## ABX (3.6e-4), alcohol (0.018)
  if (model){
    return(gg)
  } else {
    ss <- summary(aov(gg))[[1]]
    return(ss[order(ss[,'Pr(>F)'],decreasing=F),])
    # ss <- summary(gg)$coefficients
    # return(ss[order(ss[,'Pr(>|t|)'],decreasing=F),])
  }
}


############## Regressions & summary ##############
Regs <- lapply(1:10,getReg)
names(Regs) <- FactorNames

for (i in 1:10){
  write.xlsx(Regs[[i]],
             file='AG_factor_associations.xlsx',
             sheetName = FactorNames[i],
             append = i>1)
}


### This function plots the linear predictors from regression with respect to some variable, var ###
### in so doing, helps us visualize the effects of focal variable when controlling for others. ###
varplot <- function(nn,var,preds=T,...){
  y <- factorILR(nn)
  gg <- glm(y~plants+age+bmi+alcohol+sex+ABX+ibd)
  ix <- setdiff(1:length(y),gg$na.action)
  if (preds){
    plot(var[ix],predict(gg),cex=2,cex.lab=2,cex.axis=2,cex.main=2,...)
  } else {
    plot(var,y,cex=2,cex.lab=2,cex.axis=2,cex.main=2,...)
  }
}

tiff('PhyCA_associations.tiff',height=1200,width=1000)
  
  par(mfrow=c(2,2))
  abx <- ABX
  abx[abx=='Unspecified'] <- NA
  abx <- factor(abx)
  varplot(1,abx,xlab='ABX',ylab='predictor',main='Firmicutes Clade',pch=16,col=lgnd$colors[1],lwd=4)
  ag <- age
  ag[ag<22] <- NA
  varplot(2,ag,xlab='Age',ylab='predictor',main='Lachnospiraceae',pch=16,col=lgnd$colors[2],ylim=c(-10,25))
  points(ag,factorILR(2),pch=16)
  varplot(4,sex,xlab='sex',ylab='predictor',main='Ruminococcaceae',pch=16,lwd=4,col=lgnd$colors[4])
  varplot(6,ibd,xlab='ibd',ylab='predictor',main='Prevotella copri',pch=16,lwd=4,col=lgnd$colors[6])
dev.off()
