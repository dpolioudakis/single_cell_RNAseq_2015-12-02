# Make graphs to choose soft-thresholding power

print("#######################################################################")
print("Starting allen-soft-thresholding-power.R script...")
rm(list=ls())
sessionInfo()

library(WGCNA)
library(cluster)
library(flashClust)

options(stringsAsFactors=FALSE)
allowWGCNAThreads()
disableWGCNAThreads() 

exDatDF <- read.csv("../data/CellCycleNormalizedExpr.csv"
                    , header = TRUE, row.names = 1)

outPathGraph <- "../analysis/graphs/WGCNA_1_power_CCnormFPM.pdf"

plotTitle <- "Scale independence"

# Choose a set of soft-thresholding powers
powers = c(1:30)
################################################################################

# Choosing the soft-thresholding power: analysis of network topology

# Call the network topology analysis function

Make_SFT_Plots <- function(sft, outPathGraph, plotTitle) {
  
  # Plot the results:
  pdf(outPathGraph, height = 10, width = 18)
  par(mfrow = c(1, 2))
  cex1 = 0.9
  
  # Scale-free topology fit index as a function of the soft-thresholding power
  plot(sft$fitIndices[ ,1], -sign(sft$fitIndices[ ,3]) * sft$fitIndices[ ,2]
       , xlab = "Soft Threshold (power)"
       , ylab = "Scale Free Topology module Fit,signed R^2", type="n"
       , main = paste(plotTitle))
  text(sft$fitIndices[ ,1], -sign(sft$fitIndices[ ,3]) * sft$fitIndices[ ,2]
       , labels = powers, cex = cex1, col = "red")
  # this line corresponds to using an R^2 cut-off of h
  abline(h = 0.90, col = "red")
  abline(h = 0.80, col = "blue")
  abline(h = 0.70, col = "orange")
  abline(h = 0.60, col = "green")
  
  # Mean connectivity as a function of the soft-thresholding power
  plot(sft$fitIndices[ ,1], sft$fitIndices[ ,5]
       , xlab = "Soft Threshold (power)"
       , ylab = "Mean Connectivity", type = "n"
       , main = paste("Mean connectivity"))
  text(sft$fitIndices[ ,1], sft$fitIndices[ ,5], labels = powers
       , cex = cex1, col = "red")
  
  dev.off()
}

sft <- pickSoftThreshold(t(exDatDF), powerVector = powers, verbose = 5
                         , corFnc = "bicor")
Make_SFT_Plots(sft, outPathGraph, plotTitle)
