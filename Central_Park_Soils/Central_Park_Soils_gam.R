library(phylofactor)
library(plotrix)
library(viridis)
library(ggpubr)
load('pf_soil')   #PhyloFactor object (old version) from Washburne et al. 2017
Data <- pf$Data
tree <- pf$tree
X <- pf$X
Tax$taxonomy <- as.character(Tax$taxonomy)
pf.gam <- PhyloFactor(Data,tree,X,frmla=Data~s(pH)+s(C)+s(N),
                      method='gam',nfactors=4,ncores=7,choice='var')

pf.gam

## The default method='gam', with the formula above, is equivalent to the following
## call below with the customized choice.fcn and choice.fcn.dependencies inputs:
# GAM <- function(y,X,PF.output=FALSE,...){
#   dataset <- cbind('Data'=y,X)
#   gg <- mgcv::gam(Data~s(pH)+s(C)+s(N),data=dataset,...)
# 
#   if (PF.output){
#     return(gg)
#     break
#   } else {
#     output <- NULL
#     output$objective <- getStats(gg)['ExplainedVar']  ## choice='var' will pull out ExplainedVar
#     output$stopStatistics <- getStats(gg)['Pval']
#     return(output)
#   }
# }
# depends <- function() library(mgcv)
# pf.gam2 <- PhyloFactor(Data,tree,X,choice.fcn = GAM,choice.fcn.dependencies = depends,
#                       nfactors=1,ncores=7)

save(list=c('pf.gam'),file='pf_gam')



nn=1
s <- summary(pf.gam,Tax,nn)
s


### let's get a stacked prediction

B.obs <- pf.BINprojection(pf.gam,rel.abund=T,prediction=F,factor = 4)
pf.gam$glms <- pf.gam$custom.output
B.pred <- pf.BINprojection(pf.gam,rel.abund=T,prediction=T,factor=4)

Data.pred <- pf.predict(pf.gam)

pp <- pf.tree(pf.gam,layout = 'rectangular',bg.alpha=1,GroupList=bins(pf.gam$basis),
              branch.length='none',top.layer=T,top.alpha = 0.4)
pp$ggplot
ggsave('Soil_tree.png',height=6,width=2)

pH <- pf.gam$X$pH
ix <- order(pH)
pH <- sort(pH)
cols <- pp$legend$colors
tiff('Soil_bins.png',height=800,width=1400)
  par(mfrow=c(1,2))
  stackpoly(pH,t(B.obs$Data[,ix]),stack=T,cex.lab=2,cex.axis=2,cex.main=2,
            main='Observations',xlab='pH',ylab='Rel. Abundance',ylim=c(0,1),
            col = cols)
  stackpoly(pH,t(B.pred$Data[,ix]),stack=T,cex.lab=2,cex.axis=2,cex.main=2,
            main='Preditions',xlab='pH',ylab='Rel. Abundance',ylim=c(0,1),
            col=cols)
dev.off()
