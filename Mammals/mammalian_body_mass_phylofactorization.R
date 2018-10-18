library(phylofactor)
# source('R/visualization_fcns.R')
library(tidyverse)
library(plyr)
library(rotl)
library(viridis)

# Reading data ------------------------------------------------------------

panth=read.csv("~/PhyloFactor/generalized_phylofactorization/Paper/Data/pantheria.csv",header=TRUE,na.strings=-999)
taxonomy <- select(panth,'Order'=MSW05_Order,
                   'Family'=MSW05_Family,
                   'Genus'=MSW05_Genus,
                   'Species'=MSW05_Species,
                   'ID'=MSW05_Binomial)

## trim to mass only
mass=panth %>% filter(!is.na(X5.1_AdultBodyMass_g)) %>%
  select('ID'=MSW05_Binomial,'mass'=X5.1_AdultBodyMass_g) %>%
  mutate('ID'=as.character(ID))

## search with rotl
ROTL=tnrs_match_names(names=mass$ID,context_name="Animals",do_approximate_matching=F) %>%
  filter(!is.na(ott_id),!is.na(unique_name))

# bad ott_ids; ott_ids were not found in tol_induced_subtree and last 5 are singleton nodes.
bad_ott_ids <- c("513426","3611121","1055363","394961","307874",
                 "211375","238425","1080286","342715","276756",
                 "112409","417950","406006","1034143","1068214",
                 "19010","563163","791705")
ROTL <- ROTL %>% filter(!ott_id %in% bad_ott_ids)

## get tree
tree=tol_induced_subtree(ott_ids=ROTL$ott_id)

## remove ott information from the tips
tree$tip.label=strip_ott_ids(tree$tip.label)

## make tip labels match ROTL
ROTL$species=gsub(" ","_",ROTL$unique_name)
mass$ID <- gsub(' ','_',mass$ID)
taxonomy$ID <- gsub(' ','_',taxonomy$ID)

## manually change species which doesn't match:
ROTL$species=revalue(ROTL$species,
                     c("Capra_hircus_(species_in_domain_Eukaryota)"="Capra_hircus"))
ROTL <- ROTL[match(tree$tip.label,ROTL$species),] %>%
  filter(species %in% mass$ID)

mass <- mass[match(ROTL$species,mass$ID),]

rm(list=c('ROTL','panth'))
tree <- drop.tip(tree,setdiff(tree$tip.label,mass$ID))
mass <- mass[match(tree$tip.label,mass$ID),]
all.equal(tree$tip.label,mass$ID)
Z <- log(mass$mass)
names(Z) <- mass$ID


taxonomy <- taxonomy %>% filter(ID %in% mass$ID)
Tax <- data.frame('ID'=taxonomy$ID,'taxonomy'=apply(taxonomy,1,paste,collapse='; '))

# Phylofactorization ------------------------------------------------------

pf <- twoSampleFactor(Z,tree,nfactors=5,ncores=7)
saveRDS(pf,'pf_BodyMass')
# pf.W <- twoSampleFactor(Z,tree,nfactors=5,method = 'Wilcox',ncores=7) 
## example of an alternative objective function - Wilcox test
pf

nn=5
s <- summary(pf,Tax,nn)
s$taxa.split[[1]]
exp(mean(Z[pf$groups[[nn]][[1]]]))

spp <- tree$tip.label[pf$groups[[nn]][[1]]]
tree$node.label[getMRCA(tree,spp)-Ntip(tree)]

## factor 1: Perissodactyla, Cetacea, Artiodactyla     = Euungulata
## factor 2: Carnivora, Pholidota                      = Ferae
## factor 3: Chiroptera, Erinaceaomorpha, Soricomorpha = Laurasiatheria
## factor 4: Pedetidae, Anomaluridae, Dipodidae, Nesomyidae,
##           Muridae, Cricetidae, Calomyscidae, Spalacidae,  
##           Castoridae, Geomyidae, Heteromyidae       = Rodent suborders Myodonta, Anomaluromorpha, Castorimorpha
## factor 5: Cercopithecidae, Hylobatidae, Hominidae   = Simian parvorder Catarrhini
## factor 6:                                           = Pteropodidae (flying foxes)
## factor 7:                                           = Arctoidea    (bears & pinnipeds)
## factor 8:                                           = Macropodidae (kangaroos)
## factor 9:                                           = Paenungulata (elephants & manatees)
## factor 10:                                          = Mysticeti  (baleen whales)

##visualizing tree:
pf$tree$edge.length <- rep(0.1,Nedge(tree))
pf.tree(pf,top.layer=T,top.alpha = 0.1)
ggsave('mammalian_body_mass10.png',height=8,width=8)

pp=pf.tree(pf,factors = 1:5,alphas=1,layout = 'rectangular')
pp$ggplot
ggsave('mammalian_body_mass5.png',height=8,width=3)
