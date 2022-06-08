### Panzootic Project
## ME Mabry
## February 24 2022
## Following: http://blog.phytools.org/2013/04/using-makesimmap-on-set-of-trees.html

## load libraries
library(phytools)
#packageVersion("phytools") # 0.7.80
library(devtools)
library(ape)


## Set working directory
setwd("/blue/soltis/kenziemabry/Covid_ToL/")

## set args
#args = commandArgs(trailingOnly=TRUE)
#dataset <- args[1] #example : TopoFree_ND_random_100trees.RData

#load(dataset)


# ReadTrees -----------------------------------------------------------------
#load("TopoCons_FBD_MCC_SpeciesTree.RData")
#load("TopoCons_ND_MCC_SpeciesTree.RData")
#load("TopoFree_FBD_MCC_SpeciesTree.RData")
#load("TopoFree_ND_MCC_SpeciesTree.RData")


# ReadTraitData -----------------------------------------------------------
# Read in trait data
#traits <- read.csv("COVID_ToL.csv")

# rename columns
#names(traits)[names(traits) == "CovidInfection..Y...yes..N...No..U...unknown..E...experimental."] <- "Infection"

# check inputs
#unique(traits$Infection)
#unique(traits$ACE2conservation)

# convert blanks to N. Has to have a value
#traits$Infection <- sub("^$", "N", traits$Infection)
#traits$ACE2conservation <- sub("^$", "N", traits$ACE2conservation)

# check traits now
#unique(traits$Infection)
#unique(traits$ACE2conservation)

# Infection
#Infection_df <- traits[ , c("tiplabel", "Infection")]
#ACE_df <- traits[ , c("tiplabel", "ACE2conservation")]

# Match data to sampling for each tree set
#Infection_df2 <- Infection_df[Infection_df$tiplabel %in% tree2$tip.label, ]
#ACE_df2 <- ACE_df[ACE_df$tiplabel %in% tree2$tip.label, ]

# Format data for simmap
#Infection_df3 <- as.vector(Infection_df2$Infection)
#names(Infection_df3) <- Infection_df2$tiplabel

#ACE_df3 <- as.vector(ACE_df2$ACE2conservation)
#names(ACE_df3) <- ACE_df2$tiplabel


# RunSimmap ---------------------------------------------------------------

# Generate 10,000 stochastic maps
#species_simmap_ACE2 <- make.simmap(tree2, ACE_df3, nsim = 10000, model="ER")

#Sys.sleep(time = 15)

# save simmap tree
#save(species_simmap_ACE2, file = "TopoFree_ND_MCC_SpeciesTree_simmap_ACE2.RData")


# PlotSimmapTree ----------------------------------------------------------

# load tree 
#load("TopoCons_FBD_MCC_SpeciesTree_simmap_INFECTION.RData")
#load("TopoCons_ND_MCC_SpeciesTree_simmap_INFECTION.RData")
#load("TopoFree_FBD_MCC_SpeciesTree_simmap_INFECTION.RData")
#load("TopoFree_ND_MCC_SpeciesTree_simmap_INFECTION.RData")

#load("TopoCons_FBD_MCC_SpeciesTree_simmap_ACE2.RData")
#load("TopoCons_ND_MCC_SpeciesTree_simmap_ACE2.RData")
#load("TopoFree_FBD_MCC_SpeciesTree_simmap_ACE2.RData")
load("TopoFree_ND_MCC_SpeciesTree_simmap_ACE2.RData")

#load("TopoCons_FBD_MCC_SpeciesTree.RData")
#load("TopoCons_ND_MCC_SpeciesTree.RData")
#load("TopoFree_FBD_MCC_SpeciesTree.RData")
load("TopoFree_ND_MCC_SpeciesTree.RData")

# Read in trait data
traits <- read.csv("COVID_ToL.csv", row.names = 1)

# rename columns
#names(traits)[names(traits) == "CovidInfection..Y...yes..N...No..U...unknown..E...experimental."] <- "Infection"

# check inputs
#unique(traits$Infection)
unique(traits$ACE2conservation)

# convert blanks to N. Has to have a value
#traits$Infection <- sub("^$", "N", traits$Infection)
traits$ACE2conservation <- sub("^$", "N", traits$ACE2conservation)

# Infection trait
#Infection_df <- traits[ , c("tiplabel", "Infection")]
ACE_df <- traits[ , c("tiplabel", "ACE2conservation")]

# Match trait DF to tree
#Infection_df2 <- Infection_df[Infection_df$tiplabel %in% tree2$tip.label, ]
ACE_df2 <- ACE_df[ACE_df$tiplabel %in% tree2$tip.label, ]


# Put taxa in the same order as tree
#Infection_df3 <- Infection_df2[match(tree2$tip.label, Infection_df2$tiplabel), ]
ACE_df3 <- ACE_df2[match(tree2$tip.label, ACE_df2$tiplabel), ]


## Summarize the stochastic maps
obj<-summary(species_simmap_ACE2)

## specify colors
#cols <- c("E" = "#D22B2B", "N" = "#D3D3D3", "Y" = "#D22B2B") #Infection

# ACE2 order: H, L, M, N, VH, VL ("HIGH" , "LOW" , "MEDIUM", "N", "VERY HIGH", "VERY LOW")
cols <- c("H" ="#FDAE61",
          "L" = "#ABD9E9", 
          "M" = "#FFFFBF", 
          "N" = "#D3D3D3", 
          "VH" = "#D7191C", 
          "VL" = "#2C7BB6")

pdf("TopoFree_ND_MCC_SpeciesTree_simmap_ACE2.pdf", height = 60, width =60)


## Calculate offset needed for tip labels
h <- max(nodeHeights(species_simmap_ACE2[[1]]))

offset.factor <- 1.001 ## increase this for greater offset

## first plot the tree with transparent color:
plotTree(rescaleSimmap(species_simmap_ACE2[[1]], model = "depth", depth = offset.factor*h), 
         color = "transparent", type = "fan", fsize = 0.05, ftype = "i", lwd = 0.1)

par(fg="transparent")

## get dimensions for legend
space <-get("last_plot.phylo",envir=.PlotPhyloEnv)

## plot one tree, adjust color, type and size
plotTree(species_simmap_ACE2[[1]], cols, type = "fan", fsize = 0.05, ftype = "i", lwd = 0.1, add = TRUE, 
         xlim = space$x.lim, ylim = space$y.lim)

## Adjust pie node color and size
par(fg="transparent")
nodelabels(pie=obj$ace, piecol = cols, cex=0.01)
tiplabels(pie= as.matrix(obj$tips),piecol=cols,cex = 0.01)
par(fg="black")

## add Infection legend: make sure this in the in the same order as the traits
#add.simmap.legend(leg = c("No Data" , "Experimental Infection" , "Natural Infection"),
#                  colors = cols,
#                  x=0.9*par()$usr[1],
#                  y=0.9*par()$usr[4],
#                  prompt = FALSE)

## add ACE2 receptor legend: make sure this in the in the same order as the traits
add.simmap.legend(leg = c("HIGH" , "LOW" , "MEDIUM", "No Data", "VERY HIGH", "VERY LOW"),
                  colors = cols,
                  x=0.9*par()$usr[1],
                  y=0.9*par()$usr[4],
                  prompt = FALSE)


dev.off()

