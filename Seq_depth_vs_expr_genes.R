# Damon Polioudakis
# 2015-12-10
# Compare sequencing depth (total reads) to number of genes detected in each
# cell

# Inputs
#   TPM calculated by Jason Stein
#   Metadata

# Outputs
#   Correlation plot of principal compenents and covariates

rm(list=ls())
sessionInfo()

library(xlsx)
library(ggplot2)
library(reshape2)

# Inputs
lfpmDatDF <- read.csv("../data/FPM.csv")
lfpmDatDF <- lfpmDatDF[ ,-c(1, 3, 4, ncol(lfpmDatDF))]

metDatDF <- read.xlsx("../metadata/PercentageReadsMapping.xlsx", 1)
################################################################################

# Format TPM and metadata so CellIDs match and filter metadata for CellIDs in TPM
lfpmCells <- gsub("X", "Cell", colnames(lfpmDatDF))
metLFPMdatDF <- metDatDF[metDatDF$CellID %in% lfpmCells, ]

# Calculate number of detected genes (TPM >= 1)
detGenes <- apply(lfpmDatDF[ ,-1], 2, function (tPMs) sum(tPMs >= 1))
# detGenes <- apply(lfpmDatDF[ ,-1], 2, function (tPMs) sum(tPMs >= 5))

# Format data into data frame for ggplot2
nGeneNreadsDF <- data.frame(NumReads = metLFPMdatDF$NumReads, detGenes
                            , VisualQC = metLFPMdatDF$VisualQC)

# Scatter plot of sequencing depth versus number of detected genes
ggplot(nGeneNreadsDF, aes(x = NumReads, y = detGenes, color = VisualQC)) +
  geom_point() +
  ylab("Number of genes TPM > 1") +
  xlab("Total Reads") +
  ggtitle("Comparing sequencing depth to number of detected genes")
ggsave("../analysis/graphs/Sequencing_depth_vs_expressed_genes.pdf")
  
# Remove outliers, plot again, and fit linear model
# Scatter plot of sequencing depth versus number of detected genes
nGnRfiltDF <- nGeneNreadsDF[nGeneNreadsDF$detGenes < 3000, ]
ggplot(nGnRfiltDF, aes(x = NumReads, y = detGenes, color = VisualQC)) +
  geom_point() +
  ylab("Number of genes TPM > 1") +
  xlab("Total Reads") +
  stat_smooth(method = "lm") +
  ggtitle("Comparing sequencing depth to number of detected genes - outliers removed")
ggsave("../analysis/graphs/Sequencing_depth_vs_expressed_genes_RmvOutliers.pdf")