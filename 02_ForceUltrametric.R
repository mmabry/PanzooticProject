### Panzootic Project
## ME Mabry
## Febuary 24 2022
# following the phytools tutorial here: https://www.phytools.org/eqg2015/asr.html

## load libraries
library(phytools)
#packageVersion("phytools") # 0.7.80
library(devtools)
library(ape)


## Set working directory
setwd("/blue/soltis/kenziemabry/Covid_ToL/")

## Only ended up using MCC trees (below) since you end of mapping on one tree anyway
# Read100randomTrees -----------------------------------------------------------------
## Read in trees
trees <- read.tree("100RandomSubsettedDatasets/TopoFree_ND_random_100trees.tre")
trees <- read.tree("100RandomSubsettedDatasets/TopoFree_FBD_random_100trees.tre")
trees <- read.tree("100RandomSubsettedDatasets/TopoCons_ND_random_100trees.tre")
trees <- read.tree("100RandomSubsettedDatasets/TopoCons_FBD_random_100trees.tre")


# ReadTraitData -----------------------------------------------------------
## Read in trait data
traits <- read.csv("COVID_ToL.csv")

# rename columns
names(traits)[names(traits) == "CovidInfection..Y...yes..N...No..U...unknown..E...experimental."] <- "Infection"

# check inputs
unique(traits$Infection)

# convert blanks to N. Has to have a value
traits$Infection <- sub("^$", "N", traits$Infection)

# Infection
Infection_df <- traits[ , c("SciName", "Infection")]

# ForceUltrametric_100trees --------------------------------------------------------
# force ultrametric
run_utlrametric <- function(treeset, out){
  utree <- list()
  for(t in 1:100){
    ## read in each tree from the set of 100 random trees
    tree <- treeset[[t]]
    
    # Match data to sampling for each tree set
    Infection_df2 <- Infection_df[Infection_df$SciName %in% tree$tip.label, ]
    
    # Format data for simmap
    Infection_df3 <- as.vector(Infection_df2$Infection)
    names(Infection_df3) <- Infection_df2$SciName
    
    # Drop backbone tip labels
    tips_to_drop <- tree$tip.label[!(tree$tip.label %in% names(Infection_df3))]
    
    # print(tips_to_drop)
    tree2 <- drop.tip(tree, tips_to_drop)

    # Make sure tree ultrametric
    utree[[t]] <- force.ultrametric(tree2)
    
  }
  
  save(utree, file = out)
  
}


# runFunction -------------------------------------------------------------
run_utlrametric(trees, outfile)


# MCC_trees ---------------------------------------------------------------
tree <- read.nexus("MCCTrees/MamPhy_fullPosterior_BDvr_Completed_5911sp_topoCons_FBDasZhouEtAl_MCC_v2_target.tre")
tree <- read.nexus("MCCTrees/MamPhy_fullPosterior_BDvr_Completed_5911sp_topoCons_NDexp_MCC_v2_target.tre")
tree <- read.nexus("MCCTrees/MamPhy_fullPosterior_BDvr_DNAonly_4098sp_topoFree_FBDasZhouEtAl_MCC_v2_target.tre")
tree <- read.nexus("MCCTrees/MamPhy_fullPosterior_BDvr_DNAonly_4098sp_topoFree_NDexp_MCC_v2_target.tre")

# TraitData ---------------------------------------------------------------
## Read in trait data
traits <- read.csv("COVID_ToL.csv")

# rename columns
names(traits)[names(traits) == "CovidInfection..Y...yes..N...No..U...unknown..E...experimental."] <- "Infection"

# check inputs
unique(traits$Infection)

# convert blanks to N. Has to have a value
traits$Infection <- sub("^$", "N", traits$Infection)

# Infection
Infection_df <- traits[ , c("tiplabel", "Infection")]

# Match data to sampling for each tree set
Infection_df2 <- Infection_df[Infection_df$tiplabel %in% tree$tip.label, ]

# Format data
Infection_df3 <- as.vector(Infection_df2$Infection)
names(Infection_df3) <- Infection_df2$tiplabel

# DropBackboneTips --------------------------------------------------------
# Drop backbone tip labels
tips_to_drop <- tree$tip.label[!(tree$tip.label %in% names(Infection_df3))]

# print(tips_to_drop)
tree2 <- drop.tip(tree, tips_to_drop)

# ForceMCCtreeUltrametric -------------------------------------------------
# Make sure tree ultrametric
tree2 <- force.ultrametric(tree2)

# save MCC species tree which is now ultrametric and has backbone species removed
save(tree2, file = "TopoFree_ND_MCC_SpeciesTree.RData")
