###04_Species level analysis using ACE
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


# Read Trees --------------------------------------------------------------

#load("TopoCons_FBD_MCC_SpeciesTree.RData")
#load("TopoCons_ND_MCC_SpeciesTree.RData")
#load("TopoFree_FBD_MCC_SpeciesTree.RData")
#load("TopoFree_ND_MCC_SpeciesTree.RData")


# Set Up Trait Data -------------------------------------------------------

# Read in trait data
traits <- read.csv("COVID_ToL.csv")

# rename columns
names(traits)[names(traits) == "CovidInfection..Y...yes..N...No..U...unknown..E...experimental."] <- "Infection"


# check inputs
unique(traits$Infection)
#unique(traits$ACE2conservation)

# convert blanks to N. Has to have a value
traits$Infection <- sub("^$", "N", traits$Infection)
#traits$ACE2conservation <- sub("^$", "N", traits$ACE2conservation)

unique(traits$Infection)
#unique(traits$ACE2conservation)

# Infection
Infection_df <- traits[ , c("tiplabel", "Infection")]
#ACE_df <- traits[ , c("tiplabel", "ACE2conservation")]

# Match data to sampling for each tree set
Infection_df2 <- Infection_df[Infection_df$tiplabel %in% tree2$tip.label, ]
#ACE_df2 <- ACE_df[ACE_df$tiplabel %in% tree2$tip.label, ]

# Format data
Infection_df3 <- as.vector(Infection_df2$Infection)
names(Infection_df3) <- Infection_df2$tiplabel

#ACE_df3 <- as.vector(ACE_df2$ACE2conservation)
#names(ACE_df3) <- ACE_df2$tiplabel


# Run and saveACE -----------------------------------------------------------------

species_ace_INFECTION <- ace(Infection_df3, tree2, model="ER",type="discrete")


save(species_ace_INFECTION, file = "TopoFree_ND_MCC_SpeciesTree_ace_INFECTION.RData")



# Plot ACE ----------------------------------------------------------------

##load corresponding tree
#load("TopoCons_FBD_MCC_SpeciesTree.RData")
#load("TopoCons_ND_MCC_SpeciesTree.RData")
#load("TopoFree_FBD_MCC_SpeciesTree.RData")
load("TopoFree_ND_MCC_SpeciesTree.RData")

## Load trait mapped tree
#load("TopoCons_FBD_MCC_SpeciesTree_ace_ACE2.RData")
#load("TopoCons_ND_MCC_SpeciesTree_ace_ACE2.RData")
#load("TopoFree_FBD_MCC_SpeciesTree_ace_ACE2.RData")
load("TopoFree_ND_MCC_SpeciesTree_ace_ACE2.RData")

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

## specify colors
#cols <- c("E" = "#D22B2B", "N" = "#D3D3D3", "Y" = "#D22B2B") #Infection

# ACE2 order: H, L, M, N, VH, VL ("HIGH" , "LOW" , "MEDIUM", "N", "VERY HIGH", "VERY LOW")
cols <- c("H" ="#FDAE61",
          "L" = "#ABD9E9", 
          "M" = "#FFFFBF", 
          "N" = "#D3D3D3", 
          "VH" = "#D7191C", 
          "VL" = "#2C7BB6")


## Plot ACE tree
pdf("TopoFree_ND_MCC_SpeciesTree_ace_ACE2.pdf", height = 60, width =60)


## Calculate offset needed for tip labels
h <- max(nodeHeights(tree2))

offset.factor <- 1.001 ## increase this for greater offset

## first plot the tree with transparent color:
plotTree(rescaleSimmap(tree2, model = "depth", depth = offset.factor*h), 
         color = "transparent", type = "fan", fsize = 0.05, ftype = "i", lwd = 0.1)

par(fg="transparent")

## get dimensions for legend
space <-get("last_plot.phylo",envir=.PlotPhyloEnv)

## plot one tree, adjust color, type and size
plotTree(tree2, cols, type = "fan", fsize = 0.05, ftype = "i", lwd = 0.1, add = TRUE, 
         xlim = space$x.lim, ylim = space$y.lim)


## Adjust pie node color and size
par(fg="transparent")

nodelabels(pie=species_ace_ACE2$lik.anc,piecol=cols,cex=0.01)

#Infection <- as.vector(Infection_df3$Infection)
#names(Infection) <- Infection_df3$SciName

ACE <- as.vector(ACE_df3$ACE2conservation)
names(ACE) <- ACE_df3$tiplabel

#tiplabels(pie=to.matrix(Infection,sort(unique(Infection))),piecol=cols,cex=0.01)
tiplabels(pie=to.matrix(ACE,sort(unique(ACE))),piecol=cols,cex=0.01)


par(fg="black")

## add  legend
add.simmap.legend(leg = c("HIGH" , "LOW" , "MEDIUM", "N", "VERY HIGH", "VERY LOW"),
                  colors = cols,
                  x=0.9*par()$usr[1],
                  y=0.9*par()$usr[4],
                  prompt = FALSE)


dev.off()
