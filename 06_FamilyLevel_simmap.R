### Panzootic Project
## ME Mabry
## February 24 2022

## load libraries
library(phytools)
#packageVersion("phytools") # 0.7.80
library(devtools)
library(ape)


## Set working directory
setwd("/blue/soltis/kenziemabry/Covid_ToL/")


# Simmap ------------------------------------------------------------------

## Load Family tree R Data files
load("TopoFree_ND_MCC_Famtree_V2.RData")

# Read in trait data
FamTraits <- read.csv("COVID_ToL_Infection.csv")
#FamTraits <- read.csv("COVID_ToL_ACE2.csv")

# check inputs
unique(FamTraits$CovidInfection)
#unique(FamTraits$ACE2conservation)

# Infection
Infection_df <- FamTraits[ , c("Family", "CovidInfection")]
#ACE_df <- FamTraits[ , c("Family", "ACE2conservation")]

# Match data to sampling for each tree set
Infection_df2 <- Infection_df[Infection_df$Family %in% famTree2$tip.label, ]
#ACE_df2 <- ACE_df[ACE_df$Family %in% famTree2$tip.label, ]

# Format data for simmap
Infection_df3 <- as.vector(Infection_df2$CovidInfection)
names(Infection_df3) <- Infection_df2$Family

#ACE_df3 <- as.vector(ACE_df2$ACE2conservation)
#names(ACE_df3) <- ACE_df2$Family

# Generate 10,000 stochastic maps for each tree
simmapTrees_MCC_Fam_Infection <- make.simmap(famTree2, Infection_df3, nsim = 10000, model="ER")



save(simmapTrees_MCC_Fam_Infection, file = "TopoCons_ND_MCC_Famtree_V2_simmap_INFECTION.Rdata")


# plot tree ---------------------------------------------------------------

## load Family MCC tree files
load("TopoCons_FBD_MCC_Famtree_V2.RData")

## Load trait mapped tree
load("TopoCons_FBD_MCC_Famtree_V2_simmap_ACE2.Rdata")


## Read in trait data
#traits <- read.csv("COVID_ToL_Infection.csv", row.names = 1)
traits <- read.csv("COVID_ToL_ACE2.csv", row.names = 1)

## check inputs
#unique(traits$CovidInfection)
unique(traits$ACE2conservation)

## traits
#Infection_df <- traits[ , c("Family", "CovidInfection")]
ACE2_df <- traits[ , c("Family", "ACE2conservation")]

## Match trait to tree
#Infection_df2 <- Infection_df[Infection_df$Family %in% famTree2$tip.label, ]
ACE2_df2 <- ACE2_df[ACE2_df$Family %in% famTree2$tip.label, ]

##put taxa in the same order as tree
#Infection_df3 <- Infection_df2[match(famTree2$tip.label, Infection_df2$Family), ]
ACE2_df3 <- ACE2_df2[match(famTree2$tip.label, ACE2_df2$Family), ]

## Summarize the stochastic maps
obj<-summary(simmapTrees_MCC_Fam_Infection)


## specify colors
# Infection order: E, N, Y,  Y_E
#cols <- c("E" = "#FDAE61", "N" = "#D3D3D3", "Y" = "#D7191C", "Y_E" = "#E85C3B")

# ACE2 order: H, H_VH, L, L_M, L_M_H, M, N, VH, VL, VL_L, VL_M
cols <- c("H" ="#FDAE61",
          "H_VH" = "#E85C3B", 
          "L" = "#ABD9E9", 
          "L_M" = "#cde8d8",
          "L_M_H" ="#CCA296", 
          "M" = "#FFFFBF", 
          "N" = "#D3D3D3", 
          "VH" = "#D7191C", 
          "VL" = "#2C7BB6", 
          "VL_L" = "#6CAAD0",
          "VL_M" = "#8BB6BA")

#cols in species tree : https://loading.io/color/feature/RdYlBu-5/
#cols <- c("HIGH" = "#FDAE61", 
#          "LOW" = "#ABD9E9",
#          "MEDIUM" = "#FFFFBF",
#          "N" = "#D3D3D3",
#          "VERY HIGH" = "#D7191C",
#          "VERY LOW" = "#2C7BB6") #ACE2


## set up pdf
pdf("TopoCons_FBD_MCC_Famtree_V2_simmap_ACE2_V2.pdf", height = 60, width =60)

## Calculate offset needed for tip labels
h <- max(nodeHeights(simmapTrees_MCC_Fam_Infection[[1]]))

offset.factor <- 1.01 ## increase this for greater offset

## first plot the tree with transparent color:
plotTree(rescaleSimmap(simmapTrees_MCC_Fam_Infection[[1]], model = "depth", depth = offset.factor*h), 
         color = "transparent", type = "fan", fsize = 2, ftype = "i", lwd = 0.1)

par(fg="transparent")

## get dimensions for legend
space <-get("last_plot.phylo",envir=.PlotPhyloEnv)

## plot one tree, adjust color, type and size
plotTree(ladderize.simmap(simmapTrees_MCC_Fam_Infection[[1]]), cols, type = "fan", fsize = 2, ftype = "i", lwd = 0.1, add = TRUE, 
         xlim = space$x.lim, ylim = space$y.lim)

## Adjust pie node color and size
par(fg="transparent")
nodelabels(pie=obj$ace, piecol = cols, cex=0.2)
tiplabels(pie= as.matrix(obj$tips),piecol=cols, cex = 0.2)
par(fg="black")

## add  legend
#add.simmap.legend(leg = c("Experimental", "No Data", "Natural Infection", "Natural and Experimental Infection"),
#                  colors = cols,
#                  x=0.9*par()$usr[1],
#                  y=0.9*par()$usr[4],
#                  prompt = FALSE)

add.simmap.legend(leg = c("High", "High_Very High", "Low", "Low_Medium", "Low_Medium_High", "Medium", "No Data", "Very High", "Very Low", "Very Low_Low", "Very Low_Medium"),
                  colors = cols,
                  x=0.9*par()$usr[1],
                  y=0.9*par()$usr[4],
                  prompt = FALSE)



dev.off()



