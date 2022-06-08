#### Panzootic Project *note- not used in final analyses. Used MCC tree due to mapping on a single tree
### ME Mabry
## Febuary 24 2022


## load libraries
library(phytools)
#packageVersion("phytools") # 0.7.80
library(devtools)
library(ape)

## Set working directory
setwd("/blue/soltis/kenziemabry/Covid_ToL/")


# TreeInfo ----------------------------------------------------------------

## A bit about the trees: from https://doi.org/10.1371/journal.pbio.3000494 and downloaded from https://data.vertlife.org/
# Patch clades were estimated with 
# (1) DNA sampling only and no topology constraints (4,098 species, called TopoFree) or 
# (2) DNA sampling plus taxonomic constraints from the global ML tree to add the remaining unsampled species (5,911 species total, called TopoCons).
# Two types of backbones were constructed: 
# (1) ND (node-dated), using 17 fossil calibrations and one root constraint from Benton et al. 2015, see section 6 of Upham et al. 2019 for list of fossil calibrations)
# (2) tip-dated FBD (fossilized birth-death), using the morphological data set of Zhou and et al. 2013 trimmed to 76 fossil (mostly Mesozoic fossils, 66 - 252 Ma) and 22 extant taxa.

## read in concatenated .tre file of 10,000 credible trees
TopoFree_ND.tree <- read.tree("DNAonly_4098sp_topoFree_NDexp/MamPhy_BDvr_DNAonly_4098sp_topoFree_ND_10000.tre")
TopoFree_FBD.tree <- read.tree("DNAonly_4098sp_topoFree_FBDasZhouEtAl/MamPhy_BDvr_DNAonly_4098sp_topoFree_FBD_10000.tre")
TopoCons_ND.tree <- read.tree("Completed_5911sp_topoCons_NDexp/MamPhy_BDvr_Completed_5911sp_topoCons_ND_10000.tre")
TopoCons_FBD.tree <- read.tree("Completed_5911sp_topoCons_FBDasZhouEtAl/MamPhy_BDvr_Completed_5911sp_topoCons_FBD_10000.tre")


# 100RandomTrees ----------------------------------------------------------
## take a random 100 of these trees (all that is needed and makes it possible to work with this dataset)
# Approaches that apply Rubin's rules to address missing data in traits and phylogenetic sampling are particularly promising, suggesting that sampling 50-100 trees is sufficient to meaningfully capture parameter uncertainty (Nakagawa et al. 2019).
# method: http://blog.phytools.org/2013/04/picking-tree-or-set-of-trees-at-random.html
TopoFree_ND_random_100trees <- sample(TopoFree_ND.tree,size=100)
TopoFree_FBD_random_100trees <- sample(TopoFree_FBD.tree,size=100)
TopoCons_ND_random_100trees <- sample(TopoCons_ND.tree,size=100)
TopoCons_FBD_random_100trees <- sample(TopoCons_FBD.tree,size=100)


# WriteTrees --------------------------------------------------------------
## Write out these trees for use
write.tree(TopoFree_ND_random_100trees, "100RandomSubsettedDatasets/TopoFree_ND_random_100trees.tre")
write.tree(TopoFree_FBD_random_100trees, "100RandomSubsettedDatasets/TopoFree_FBD_random_100trees.tre")
write.tree(TopoCons_ND_random_100trees, "100RandomSubsettedDatasets/TopoCons_ND_random_100trees.tre")
write.tree(TopoCons_FBD_random_100trees, "100RandomSubsettedDatasets/TopoCons_FBD_random_100trees.tre")
