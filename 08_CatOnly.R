### Panzootic Project
## ME Mabry
## Febuary 24 2022
## Useful Tutorials:
# http://blog.phytools.org/2013/04/using-makesimmap-on-set-of-trees.html
# https://www.phytools.org/eqg2015/asr.html
# http://blog.phytools.org/2015/11/describesimmap-when-summarizing-results.html

## load libraries
library(phytools)
#packageVersion("phytools") # 0.7.80
library(devtools)
library(ape)
library(RColorBrewer)
library(geiger)

## Set working directory
setwd("/blue/soltis/kenziemabry/Covid_ToL/")


# load MCC trees ----------------------------------------------------------

## load MCC  tree files 
tree <- read.nexus("MCCTrees/MamPhy_fullPosterior_BDvr_DNAonly_4098sp_topoFree_NDexp_MCC_v2_target.tre")


# extract Felidae clade -----------------------------------------------------------

## Read in trait data to get Family names
traits <- read.csv("COVID_ToL.csv")

## new names
names_df <- traits[ , c("tiplabel", "Family")]

cat_df <- names_df[names_df$Family == "FELIDAE", ]

## Match to tree
cat_df2 <- cat_df[cat_df$tiplabel %in% tree$tip.label, ]


## make vector of species to keep 
cats_to_keep <- as.vector(cat_df2$tiplabel)

## keep only cats
CatTree <- keep.tip(tree, cats_to_keep)
plot(CatTree)

## save new tree
save(CatTree, file = "TopoFree_ND_CatTree.RData")


# run simmap --------------------------------------------------------------

load("TopoCons_FBD_CatTree.RData")

### Get trait data
##rename columns
names(traits)[names(traits) == "CovidInfection..Y...yes..N...No..U...unknown..E...experimental."] <- "Infection"

## get characters
traits$Infection <- sub("^$", "N", traits$Infection)

Infection_df <- traits[ , c("tiplabel", "Infection")]

unique(traits$Infection)

# Match data to sampling for each tree set
Infection_df2 <- Infection_df[Infection_df$tiplabel %in% CatTree$tip.label, ]


# Format data for simmap
Infection_df3 <- as.vector(Infection_df2$Infection)
names(Infection_df3) <- Infection_df2$tiplabel


# Generate 10,000 stochastic maps
simmapTrees_MCC_CAT_Infection <- make.simmap(CatTree, Infection_df3, nsim = 10000, model="ER")

save(simmapTrees_MCC_CAT_Infection, file = "TopoFree_FBD_MCC_CatTree_simmap_Infection.Rdata")



# run ace -----------------------------------------------------------------

## ACE method
aceTrees_MCC_CAT_Infection <- ace(Infection_df3, CatTree, model="ER", type="discrete")

save(aceTrees_MCC_CAT_Infection, file = "TopoCons_FBD_MCC_CatTree_ace_Infection.Rdata")


# plot simmap -------------------------------------------------------------

## Plot simmap tree
load("TopoFree_FBD_CatTree.RData")
load("TopoFree_FBD_MCC_CatTree_simmap_Infection.Rdata")

catTraits <- read.csv("COVID_ToL_CAT.csv")

names(catTraits)[names(catTraits) == "CovidInfection..Y...yes..N...No..U...unknown..E...experimental."] <- "Infection"

## get characters
catTraits$Infection <- sub("^$", "N", catTraits$Infection)

Infection_df <- catTraits[ , c("tiplabel", "Infection")]

Infection_df2 <- Infection_df[Infection_df$tiplabel %in% CatTree$tip.label, ]

##put taxa in the same order as tree
Infection_df3 <- Infection_df2[match(CatTree$tip.label, Infection_df2$tiplabel), ]

## Summarize the stochastic maps
obj<-summary(simmapTrees_MCC_CAT_Infection)

## specify colors
# Infection order: N, Y
cols <- c("N" = "#D3D3D3", "Y" = "#D22B2B")

## set up pdf
pdf("TopoFree_FBD_MCC_CatTree_simmap_INFECTION_V2.pdf", height = 30, width =20)

## Calculate offset needed for tip labels
#h <- max(nodeHeights(simmapTrees_MCC_CAT_Infection[[1]]))

#offset.factor <- 1.01 ## increase this for greater offset

## first plot the tree with transparent color:
#plotTree(rescale(simmapTrees_MCC_CAT_Infection[[1]], model = "depth", depth = offset.factor*h), 
#         color = "transparent", fsize = 1.5, ftype = "i", lwd = 0.5)

#par(fg="transparent")

## get dimensions for legend
#space <-get("last_plot.phylo",envir=.PlotPhyloEnv)

## plot one tree, adjust color, type and size
#plotTree(simmapTrees_MCC_CAT_Infection[[1]], cols, fsize = 1.5, ftype = "i", lwd = 0.5, add = TRUE, 
#         xlim = space$x.lim, ylim = space$y.lim)
plotTree(simmapTrees_MCC_CAT_Infection[[1]], cols, fsize = 2, ftype = "i", lwd = 1)

## Adjust pie node color and size
par(fg="transparent")
nodelabels(pie=obj$ace, piecol = cols, cex=0.5)
tiplabels(pie= as.matrix(obj$tips),piecol=cols, cex = 0.5)
par(fg="black")

## add  legend
add.simmap.legend(leg = c("No Data", "Natural Infection"),
                  colors = cols,
                  #x=0.9*par()$usr[1],
                  #y=0.9*par()$usr[4],
                  prompt = FALSE)


dev.off()


# plot ace ----------------------------------------------------------------

## plot ace tree
load("TopoFree_FBD_CatTree.RData")
load("TopoFree_FBD_MCC_CatTree_ace_Infection.Rdata")

catTraits <- read.csv("COVID_ToL_CAT.csv")

names(catTraits)[names(catTraits) == "CovidInfection..Y...yes..N...No..U...unknown..E...experimental."] <- "Infection"

## get characters
catTraits$Infection <- sub("^$", "N", catTraits$Infection)

Infection_df <- catTraits[ , c("tiplabel", "Infection")]

Infection_df2 <- Infection_df[Infection_df$tiplabel %in% CatTree$tip.label, ]

##put taxa in the same order as tree
Infection_df3 <- Infection_df2[match(CatTree$tip.label, Infection_df2$tiplabel), ]

## specify colors
# Infection order: N, Y
cols <- c("N" = "#D3D3D3", "Y" = "#D22B2B")

## set up pdf
pdf("TopoFree_FBD_MCC_CatTree_ace_INFECTION_V5.pdf", height = 30, width =20)

## Calculate offset needed for tip labels
#h <- max(nodeHeights(CatTree))

#offset.factor <- 1.0001 ## increase this for greater offset

## first plot the tree with transparent color:
#plotTree(rescale(CatTree, model = "depth", depth = offset.factor*h), 
#         color = "transparent", fsize = 1.5, ftype = "i", lwd = 0.5)

#par(fg="transparent")

## get dimensions for legend
#space <-get("last_plot.phylo",envir=.PlotPhyloEnv)

## plot one tree, adjust color, type and size
#plotTree(CatTree, cols, fsize = 1.5, ftype = "i", lwd = 0.5, add = TRUE, 
#         xlim = space$x.lim, ylim = space$y.lim)

plotTree(CatTree, cols, fsize = 2, ftype = "i", lwd = 1)


## Adjust pie node color and size
par(fg="transparent")
nodelabels(pie=aceTrees_MCC_CAT_Infection$lik.anc,piecol=cols,cex=0.5)

Infection <- as.vector(Infection_df3$Infection)
names(Infection) <- Infection_df3$tiplabel

tiplabels(pie=to.matrix(Infection,sort(unique(Infection))),piecol=cols,cex=0.5)

par(fg="black")

## add  legend
add.simmap.legend(leg = c("No Data", "Natural Infection"),
                  colors = cols,
                  #x=0.9*par()$usr[1],
                  #y=0.9*par()$usr[4],
                  prompt = FALSE)

dev.off()
