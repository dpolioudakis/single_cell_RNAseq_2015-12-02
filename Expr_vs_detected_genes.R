# Damon Polioudakis
# 2015-12-12
# Examine relationship between Median Expression and Number of detected genes

# Inputs
#   HTseq counts expression
#   Metadata
#   TPM generated by Jason Stein

rm(list=ls())
sessionInfo()

library(xlsx)
library(ggplot2)
library(reshape2)

# Inputs
# TPM from Jason Stein
lfpmDatDF <- read.csv("../data/FPM.csv")
lfpmDatDF <- lfpmDatDF[ ,-c(1, 3, 4, ncol(lfpmDatDF))]
# Gene expression from HTSeq
exDatDF <- read.csv("../SxaQSEQsXap089L2_HTSC/Exprs_HTSCexon.csv", row.names = 1)
# Metadata
metDatDF <- read.xlsx("../metadata/PercentageReadsMapping.xlsx", 1)
################################################################################

# Plots compareing Median Expression versus Number of detected genes

# Caculate TPMs myself

# Edit CellIDs to be same format
colnames(exDatDF) <- gsub("_.*$", "", colnames(exDatDF))
sortMetDF <- metDatDF[order(metDatDF$CellID), ]
tpmDF <- t(t(exDatDF) / (sortMetDF$NumReads / 1000000))
# Filter genes to only those with TPM > 5 in at least 10% of samples
keep <- apply(data.frame(tpmDF > 5), 1, function(expd) sum(expd) >= 0.1*ncol(tpmDF))
filtTPMdF <- (tpmDF[keep, ])

# Calculate number of detected genes (TPM >= 1)
detGenes <- apply(filtTPMdF, 2, function (tPMs) sum(tPMs >= 1))
# detGenes <- apply(lfpmDatDF[ ,-1], 2, function (tPMs) sum(tPMs >= 5))

# Calculate median expression for each cell
medExpr <- apply(filtTPMdF, 2, function (tPMs) median(tPMs[tPMs >= 1]))

# Format data into data frame for ggplot2
nGeneMedExprDF <- data.frame(MedianExpr = medExpr, detGenes
                             , VisualQC = sortMetDF$VisualQC)

# Remove outliers, plot again, and fit loess model
nGnRfiltDF <- nGeneMedExprDF[nGeneMedExprDF$detGenes < 3000, ]
nGnRfiltDF <- nGnRfiltDF[nGnRfiltDF$MedianExpr < 3000, ]
# Scatter plot of sequencing depth versus number of detected genes
ggplot(nGnRfiltDF, aes(x = detGenes, y = MedianExpr, color = VisualQC)) +
  geom_point() +
  stat_smooth(data = nGnRfiltDF[nGnRfiltDF$VisualQC == "A", ], method = "loess") +
  stat_smooth(method = "lm") +
  ylab("Median Expression") +
  xlab("Number of detected genes ( > 1 TPM)") +
  ggtitle("Comparing median expression to number of detected genes for each cell
          TPMs calculated by Damon Polioudakis")


# Using TPMs from Jason Stein

# Format TPM and metadata so CellIDs match and filter metadata for CellIDs in TPM
lfpmCells <- gsub("X", "Cell", colnames(lfpmDatDF))
metLFPMdatDF <- metDatDF[metDatDF$CellID %in% lfpmCells, ]

# Calculate number of detected genes (TPM >= 1)
detGenes <- apply(lfpmDatDF[ ,-1], 2, function (tPMs) sum(tPMs >= 1))
# detGenes <- apply(lfpmDatDF[ ,-1], 2, function (tPMs) sum(tPMs >= 5))

medExpr <- apply(lfpmDatDF[ ,-1], 2, function (tPMs) median(tPMs[tPMs >= 1]))

# Format data into data frame for ggplot2
nGeneMedExprDF <- data.frame(MedianExpr = medExpr, detGenes
                             , VisualQC = metLFPMdatDF$VisualQC)

# Remove outliers, plot again, and fit loess model
nGnRfiltDF <- nGeneMedExprDF[nGeneMedExprDF$detGenes < 3000, ]
# Scatter plot of sequencing depth versus number of detected genes
ggplot(nGnRfiltDF, aes(x = detGenes, y = MedianExpr, color = VisualQC)) +
  geom_point() +
  stat_smooth(data = nGnRfiltDF[nGnRfiltDF$VisualQC == "A", ], method = "loess") +
  ylab("Median Expression (TPM)") +
  xlab("Number of detected genes ( > 1 TPM)") +
  ggtitle("Comparing median expression to number of detected genes for each cell
          TPMs from Jason Stein")
ggsave("../analysis/graphs/Expr_vs_detected_genes_TPM.pdf")

# Scale TPMs from Jason Stein by proportion of detected genes for each cell

# Proportion of detected genes for each cell
propDetG <- detGenes / nrow(lfpmDatDF)
nDetNrmDF <- data.frame(apply(lfpmDatDF[ ,-1], 2
                      , function (tPMs) tPMs * (propDetG / (median(propDetG)))))
medExpr <- apply(nDetNrmDF, 2, function (tPMs) median(tPMs[tPMs >= 1]))

# Format data into data frame for ggplot2
nGeneMedExprDF <- data.frame(MedianExpr = medExpr, detGenes
                             , VisualQC = metLFPMdatDF$VisualQC)

# Remove outliers, plot again, and fit loess model
nGnRfiltDF <- nGeneMedExprDF[nGeneMedExprDF$detGenes < 3000, ]
# Scatter plot of sequencing depth versus number of detected genes
ggplot(nGnRfiltDF, aes(x = detGenes, y = MedianExpr, color = VisualQC)) +
  geom_point() +
  stat_smooth(data = nGnRfiltDF[nGnRfiltDF$VisualQC == "A", ], method = "loess") +
  ylab("Median Expression (TPM scaled)") +
  xlab("Number of detected genes ( > 1 TPM)") +
  ggtitle("Comparing median expression to number of detected genes for each cell
          TPMs from Jason Stein scaled by
          (Proportion of detected genes for a cell / median proportion of detected genes for all cells)")
ggsave("../analysis/graphs/Expr_vs_detected_genes_TPMscaled.pdf")

