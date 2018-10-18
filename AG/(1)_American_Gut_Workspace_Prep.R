### American Gut phylofactorizaztion
library(utils)
library(phyloseq)
library(ape)
library(biomformat)
library(magrittr)
setwd('AG/')
# download.file(url='ftp://ftp.microbio.me/AmericanGut/ag-July-29-2016/03-otus/100nt/gg-13_8-97-percent/97_otus.tree',destfile='97_otus.tree')
# download.file(url='ftp://ftp.microbio.me/AmericanGut/ag-July-29-2016/03-otus/100nt/gg-13_8-97-percent/otu_table.biom',destfile='otu_table.biom',mode='wb')
# download.file(url='ftp://ftp.microbio.me/AmericanGut/ag-July-29-2016/04-meta/ag-cleaned.txt',destfile = 'ag-cleaned.txt')


############### OTU Table ##################
############################################
source('import_biom2.R')
data <- read_hdf5_biom('otu_table.biom') %>% import_biom2(biom) # now a phyloseq object
gc()

OTUTable <- data@otu_table@.Data
otus <- rownames(OTUTable)
samples <- colnames(OTUTable)


############### Taxonomy ##################
############################################
taxonomy <- data@tax_table@.Data
Taxonomy <- matrix(NA,ncol=2,nrow=nrow(taxonomy))
Taxonomy[,1]<-rownames(taxonomy)
Taxonomy[,2]<- apply(taxonomy,MARGIN=1,FUN=function(b) paste(b,collapse="; "))
colnames(Taxonomy) <- c('OTU.IDs','taxonomy')
rm('taxonomy')
rm('data')
gc()

############### tree ##################
############################################
tree <- read_tree_greengenes('97_otus.tree')
# there were 13 warnings: In go.down() : NAs introduced by coercion

############### meta-data ##################
############################################
X <- read.table('ag-cleaned.txt',header=T,sep='\t',quote = "",comment.char = '')
colnames(X)[1] <- 'SampleID'

############### cross-checking ############
########### and eliminating NA samples ####

all(otus %in% tree$tip.label)
#TRUE - all the otus in our table are in 


all(samples %in% X$SampleID)
#FALSE
ix <- which(samples %in% X$SampleID)
samples <- samples[ix]
OTUTable <- OTUTable[,samples]
OTUTable <- OTUTable[which(rowSums(OTUTable)>0),]
otus <- rownames(OTUTable)

X <- X[which(X$SampleID %in% colnames(OTUTable)),]
ncol(OTUTable)==nrow(X)

Taxonomy <- Taxonomy[which(Taxonomy[,1] %in% otus),]

length(tree$tip.label)-length(setdiff(tree$tip.label,otus))
tree <- drop.tip(tree,setdiff(tree$tip.label,otus))


all(rownames(OTUTable) %in% tree$tip.label)
#TRUE
all(colnames(OTUTable) %in% X$SampleID)
#TRUE
all(rownames(OTUTable) %in% Taxonomy[,1])
#TRUE

save(list=c('tree','Taxonomy','X','OTUTable'),file='AmericanGut_Workspace')
