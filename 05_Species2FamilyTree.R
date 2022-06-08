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


# ReadSpeciesMCCtrees -----------------------------------------------------

## load MCC species tree files 
#load("TopoFree_ND_MCC_SpeciesTree.RData")
#load("TopoFree_FBD_MCC_SpeciesTree.RData")
#load("TopoCons_ND_MCC_SpeciesTree.RData")
load("TopoCons_FBD_MCC_SpeciesTree.RData")

length(tree2$tip.label)


# ReadTraitData -----------------------------------------------------------

# Read in trait data to get Family names
traits <- read.csv("COVID_ToL.csv")

# new names
names_df <- traits[ , c("tiplabel", "Family")]


# RenameTipsFunction ------------------------------------------------------

## Function for renaming tips
# new_names and old_names are columns of a dataframe with the taxon names
rename.tips.phylo <- function(tree, names) {
  tree2 <- tree
  tree2$tip.label <- names$Family[match(tree2$tip.label,names$tiplabel)]
  return(tree2)
}


# MakeFamilyTree ----------------------------------------------------------

famTree <- rename.tips.phylo(tree2, names_df)


# GetUniqueFamilys --------------------------------------------------------

# get list of unique family names in tree
FamTreeTipLabels <- unique(famTree$tip.label)


# ReduceSpeciesTree2Family ------------------------------------------------

# for each family name ...
for(i in FamTreeTipLabels){
  fam <- i
  print(fam)
  
  # get possible species names for family
  drop_tips <- names_df[names_df$Family == fam, ]$tiplabel
  print(drop_tips)
  
  # find which of these species are in your tree
  drop_tips_in_tree <- drop_tips[drop_tips %in% tree2$tip.label]
  
  # if you just have one species name it goes to the next family
  if(length(drop_tips_in_tree) == 1){
    next
  }
  
  # using the species in the tree as drop list, except the first one
  tree2 <- drop.tip(tree2, drop_tips_in_tree[2:length(drop_tips_in_tree)])
}



# ReplaceSpeciesWithFamily ------------------------------------------------

## Take trees with one species per family and replace with family name
famTree2 <- rename.tips.phylo(tree2, names_df)

## save new tree
save(famTree2, file = "TopoCons_FBD_MCC_Famtree_V2.RData")
