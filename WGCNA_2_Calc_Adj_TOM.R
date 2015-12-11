#!/usr/bin/Rscript

# Damon Polioudakis
# 2015-12-07
# WGCNA: Cluster samples and construct modules

# Workflow
  # WGCNA_1_Soft_Thresholding.R
  # WGCNA_2_Calc_Adj_TOM.R

print("#######################################################################")
print("Starting WGCNA_2_Calc_Adj_TOM.R script...")
rm(list=ls())
sessionInfo()

library(WGCNA)

options(stringsAsFactors=FALSE)
allowWGCNAThreads()
# disableWGCNAThreads()
argsL <- commandArgs(TRUE)

softPower <- as.numeric(argsL[[1]])
print(paste("Soft Power:", softPower))
# softPower <- 17

# Input
exDatDF <- read.csv("../data/CellCycleNormalizedExpr.csv", row.names = 1)
exDatDF <- t(exDatDF)

# Output paths
outTOMadj <- paste("../data/WGCNA_2_Adjacency_TOM_CCnormFPM_SP", softPower
                   , ".rda", sep = "")
print("#######################################################################")

# Biweight midcorrelation is considered to be a good alternative to Pearson
# correlation since it is more robust to outliers.

print("Starting adjacency calculation...")
adjacency <- adjacency(exDatDF, power = softPower, corFnc= "bicor"
                       , type = "signed")
print("Finished adjacency calculation...")

print("Starting TOM calculation...")
TOM <- TOMsimilarity(adjacency, TOMType = "signed")
dissTOM <- 1-TOM
print("Finished TOM calculation...")

save(adjacency, TOM, dissTOM, file = outTOMadj)

print("End of WGCNA_2_Calc_Adj_TOM.R script...")