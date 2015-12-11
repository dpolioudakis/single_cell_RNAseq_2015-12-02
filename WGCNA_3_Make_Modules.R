#!/usr/bin/Rscript

# Damon Polioudakis
# 2015-12-07
# WGCNA: Cluster samples and construct modules
# Recommended script call:
  # Rscript WGCNA_3_Make_Modules.R | tee logs/WGCNA_3_Make_Modules_$(date +%Y%m%d%H%M%S).log

# Workflow
  # WGCNA_1_Soft_Thresholding.R
  # WGCNA_2_Calc_Adj_TOM.R
  # WGCNA_3_Make_Modules.R

rm(list=ls())
print("#######################################################################")
print("Starting WGCNA_3_Make_Modules.R script...")
sessionInfo()

library(WGCNA)
library(cluster)
library(flashClust)

options(stringsAsFactors=FALSE)
allowWGCNAThreads()
disableWGCNAThreads() 

load("../data/WGCNA_2_Adjacency_TOM_CCnormFPM_SP17.rda")
exDatDF <- read.csv("../data/CellCycleNormalizedExpr.csv", row.names = 1)
exDatDF <- t(exDatDF)

outModules <- "../data/WGCNA_3_Modules_CCnormFPM_SP17.rda"
outDendModGraph <- "../analysis/graphs/WGCNA_3_DendMods_CCnormFPM_SP17.pdf"
print("#######################################################################")

# Make Modules

geneTree <- hclust(as.dist(dissTOM), method = "average")

# Module identification using hybrid tree cut:
# Function to construct modules
# Args: minimum module size, deep split
Make_modules <- function (minModSize, deepSplit) {
  print("Treecut arguments:")
  # print(c("minModSize"=minModSize,"cut.height"=cutHeightMergeME, "deepSplit"=ds))
  print(c(minModSize, deepSplit))
  tree = cutreeDynamic(dendro = geneTree, pamRespectsDendro = TRUE
                      , minClusterSize = minModSize
                      # , cutHeight = 0.967
                      , deepSplit = deepSplit, distM = as.matrix(dissTOM))
  tree <- labels2colors(tree)
  print("Table of genes per module:")
  print(table(tree))
  tree
}
# cutreeHyrid make modules function
# Make_modules <- function (minModSize, deepSplit) {
#   print("Treecut arguments:")
#   # print(c("minModSize"=minModSize,"cut.height"=cutHeightMergeME, "deepSplit"=ds))
#   print(c(minModSize, deepSplit))
#   tree = cutreeHybrid(dendro = geneTree, pamRespectsDendro = TRUE
#                       , minClusterSize = minModSize
#                       # , cutHeight = 0.967
#                       , deepSplit = deepSplit, distM = as.matrix(dissTOM))
#   print("Table of genes per module:")
#   print(table(tree$labels))
#   tree$labels
# }

# Merge modules based on ME function
# Args: Modules colors, Cut height to merge ME
Merge_modules_ME <- function (genesModuleColor, cutHeightMergeME) {
  # Call an automatic merging function
  # merged: The merged module colors
  # Cut height of 0.25, corresponds to a correlation of 0.75, to merge ME:
  merged <- mergeCloseModules(exprData = exDatDF, colors = genesModuleColor,
                              cutHeight = cutHeightMergeME)
  print("Table of genes per module after merging:")
  print(table(merged$colors))
  merged$colors
}
# For cutreehybrid
# Merge_modules_ME <- function (genesModuleColor, cutHeightMergeME) {
#   # Call an automatic merging function
#   # merged: The merged module colors
#   # Cut height of 0.25, corresponds to a correlation of 0.75, to merge ME:
#   merged <- mergeCloseModules(exprData = exDatDF, colors = genesModuleColor,
#                               cutHeight = cutHeightMergeME)
#   print("Table of genes per module after merging:")
#   print(table(merged$colors))
#   labels2colors(merged$colors)
# }

# Test different parameters for constructing and merging modules
# Define arguments to test for cutreeHybrid

minModSizes <- c(30, 50, 100)
deepSplits <- c(2, 4)
cutHeightMergeMEs <- c(0, 0.1, 0.2, 0.25)

modulesColors <- NULL
moduleParameterLabels <- NULL
for (minModSize in minModSizes) {
  for (deepSplit in deepSplits) {
    # Test multiple cutreeHybrid parameters
    module <- Make_modules(minModSize, deepSplit)
    for (cutHeightMergeME in cutHeightMergeMEs) {
      # Test ME merge cut heights
      modulesColors <- cbind(modulesColors,
                             Merge_modules_ME(module, cutHeightMergeME))
      # Make label from parameters used to make each module
      moduleParameterLabels <- c(moduleParameterLabels, paste(
        "MMS=",minModSize
        , " \nDS=",deepSplit
        , " \nMEcor=",cutHeightMergeME
      ))
    }
  }
}
print("Done making modules...")

#sizeGrWindow(25,20)
print("Making WGCNA modules dendrogram graph...")
pdf(file = outDendModGraph, height = 25, width = 20)
plotDendroAndColors(geneTree
                    , modulesColors
                    , groupLabels = moduleParameterLabels
                    , addGuide = TRUE
                    , dendroLabels = FALSE
                    , main = "Dendrogram With Different Module Cutting Parameters")
dev.off()

save(exDatDF, geneTree, modulesColors, moduleParameterLabels, file = outModules)
