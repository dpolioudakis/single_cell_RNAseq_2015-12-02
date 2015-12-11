rm(list=ls())

source("http://bioconductor.org/biocLite.R")
biocLite("Genominator")
library(Genominator)

options(stringsAsFactors=FALSE)

exDatDF <- read.csv("../SxaQSEQsXap089L2_HTSC/Exprs_HTSCexon.csv")

## Get Gencode 18 gtf file - this was cleaned by selecting only the columns containing the word "exon" and with the relevant information - chrnum, feature type, strand, start, end, and ensg and ense IDs separated by a semicolon
gtfinfo <- read.table("../source/gencode.v19.annotation.gtf",sep="\t")
keep <- gtfinfo[,3]=="exon" ## Keep only the exon level features
gtfinfo <- gtfinfo[keep,]

genexoninfo <- unlist(strsplit(gtfinfo[,9],"[;]")) ## Split the semicolon separated information
gen.col<- which(regexpr("gene_id ", genexoninfo)>0)  ## finding which has gene_id for ensembl ID
getseq= genexoninfo[gen.col]
ENSGID <- substr(getseq,9,100)
length(unique(ENSGID)) ##62069

trans.col=which(regexpr("transcript_id ", genexoninfo)>0)
transeq = genexoninfo[trans.col]
ENSEID <- substr(transeq,16,100)
length(unique(ENSEID)) #213272
gtfinfo <- cbind(gtfinfo[,c(1:8)],ENSGID,ENSEID)

gene.col=which(regexpr("gene_name ", genexoninfo)>0)
geneseq = genexoninfo[gene.col]
GENEID <- substr(geneseq,12,100)
length(unique(GENEID)) #52775


gtfinfo <- cbind(gtfinfo[,c(1:8)],ENSGID,ENSEID)
gtfinfo= gtfinfo[,-c(6,8)] ## 6 and 8 columns are blank

## Keep only one copy of each ENSEID - the gtf file records one copy for each transcript id
keep <- match(unique(ENSEID),ENSEID)
gtfinfo1 <- gtfinfo[keep,]
##gtfinfo[,1] <- substr(gtfinfo[,1],4,10) ## 672406 exons is exactly what biomaRt contains

## Recode things for the Genominator package
chrnums <- gtfinfo1[,1] ## Using as.factor to coerce chromosome names can really botch things up... beware! So go ahead and convert MT, X, and Y to numbers throughout, unless necessary for other purposes
chrnums[chrnums=="MT"] <- "20"
chrnums[chrnums=="X"] <- "21"
chrnums[chrnums=="Y"] <- "22"
# rmChR.col1=which(regexpr("HG", chrnums)>0)
# rmChR.col2= which(regexpr("GL", chrnums)>0) ## removing Non-annotated(NT) chromosomes
# rmChR.col3= which(regexpr("HS", chrnums)>0)
# rmChR.col=c(rmChR.col1,rmChR.col2,rmChR.col3)
# gtfinfo1=gtfinfo1[-rmChR.col,]
chrnums=chrnums[-rmChR.col]
gtfinfo1[,1] <- chrnums ## Check here

gtfinfo=gtfinfo1

strinfo <- gtfinfo[,6]
strinfo[strinfo=="+"] <- 1L
strinfo[strinfo=="-"] <- -1L
gtfinfo[,6] <- strinfo


geneDat1=gtfinfo[,c(1,6,4,5,7,8)] ## chr integer, strand integer (-1L,0L,1L), start integer, end integer, ensg and transcript id
geneDat1 <- data.frame(as.numeric(chrnums),as.numeric(geneDat1[,2]),as.numeric(geneDat1[,3]),as.numeric(geneDat1[,4]),geneDat1[,5],geneDat1[,6])
names(geneDat1) <- c("chr","strand","start","end","ensembl_gene_id","ensembl_exon_id")

geneDatX <- geneDat1[order(geneDat1[,1],geneDat1[,3]),]
#DP edit - remove NAs from ERCC chromosomes
geneDatX <- geneDatX[complete.cases(geneDatX),]
validAnnotation(geneDatX) ## Have genominator check if this is a valid data object
geneDatX <- makeGeneRepresentation(annoData=geneDat1,type="Ugene",gene.id = "ensembl_gene_id", transcript.id = "ensembl_exon_id",verbose=TRUE) ##should take few minutes !!!

save(geneDatX,file="../source/GenominatorUnionGeneModelsENSEMBLhg19.rda")
load(file="../source/GenominatorUnionGeneModelsENSEMBLhg19.rda")

## Now use the genominator output to calculate GC content ### Use mac laptop ###
library(Repitools)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg19)
library(GenomicRanges)
geneDat2 <- cbind(geneDatX,geneDatX[,3]-geneDatX[,2])
geneDat2 <- geneDat2[order(geneDat2[,5]),]

## Change formatting again
chrnums <- geneDat2[,"chr"]
chrnums[chrnums=="20"] <- "M" ## important as UCSC codes MT as M
chrnums[chrnums=="21"] <- "X"
chrnums[chrnums=="22"] <- "Y"
stinfo <- geneDat2[,"strand"]
stinfo[stinfo==c(-1)] <- "-"
stinfo[stinfo==c(1)] <- "+"

## Calculate GC content from hg19 using the union exon ranges
gcQuery <- GRanges(paste("chr", chrnums,sep=""),IRanges(geneDat2[,2],geneDat2[,3]),strand=stinfo) ## Convert to a genomic ranges object
gcContent <- gcContentCalc(x=gcQuery,organism=Hsapiens)

## Take a length weighted average of GC content percentages to get the GC content for the union gene model
head(geneDat2)
geneDat2 <- cbind(geneDat2,gcContent)
geneDat2 <- cbind(geneDat2,gcContent*geneDat2[,6])
unionGenes <- by(geneDat2[,6],as.factor(geneDat2[,5]),sum)
unionGC <- by(geneDat2[,8],as.factor(geneDat2[,5]),sum)
geneDat3 <- cbind(unionGenes,unionGC/unionGenes)
colnames(geneDat3) <- c("UnionGeneLength","UnionGCcontent")
ENSEMBLhg19.70UnionAnno <- geneDat3

## Save for further usage
save(ENSEMBLhg19.70UnionAnno, file="../source/ENSEMBLhg19_UnionAnno.rda")
load("../source/ENSEMBLhg19_UnionAnno.rda")

# DP edit - calculate length bias for each cell
unionGenes <- cbind(Length = unionGenes)
exLenDF <- merge(x = exDatDF, y = unionGenes, by.x = "X", by.y = "row.names" )
lenBias <- apply(exLenDF, 2
                 , function(counts) sum(as.numeric(counts) * exLenDF["Length"]) / 
                   sum(as.numeric(counts)))
gcBiasDF <- merge(x = exDatDF, y = ENSEMBLhg19.70UnionAnno, by.x = "X", by.y = "row.names" )
gcBiasDF <- gcBiasDF[complete.cases(gcBiasDF), ]
gcBias <- apply(gcBiasDF, 2
                 , function(counts) sum(as.numeric(counts) * gcBiasDF["UnionGCcontent"]) / 
                   sum(as.numeric(counts)))
save(lenBias, gcBias, file = "../analysis/tables/Gene_Length_GC_Bias.rda")







## Check difference in relationship to GC before and after normalization
preNorm <- datExpr.HTSC
postNorm <- cqn.dat
keepgenes <- intersect(rownames(preNorm),rownames(postNorm))
preNorm <- preNorm[match(keepgenes,rownames(preNorm)),]
postNorm <- postNorm[match(keepgenes,rownames(postNorm)),]
geneAnno1 <- geneAnno[match(keepgenes,rownames(geneAnno)),]

qualMat <- matrix(NA,nrow=ncol(preNorm),ncol=4)
colnames(qualMat) <- c("pre.Norm.GC.cor","pre.Norm.Length.cor","post.Norm.GC.cor","post.Norm.Length.cor")

for (i in 1:nrow(qualMat)) {
  qualMat[i,1] <- cor(preNorm[,i],geneAnno1[,2],method="spearman")
  qualMat[i,2] <- cor(preNorm[,i],geneAnno1[,1],method="spearman")
  qualMat[i,3] <- cor(postNorm[,i],geneAnno1[,2],method="spearman")
  qualMat[i,4] <- cor(postNorm[,i],geneAnno1[,1],method="spearman")
}
quantile(qualMat[,1])
quantile(qualMat[,3])
quantile(qualMat[,2])
quantile(qualMat[,4])
pdf("GCcorrelations.pdf",width=8,height=8)
par(mfrow=c(2,2))
hist(qualMat[,1],main=colnames(qualMat)[1],xlim=c(-0.15,0.1),ylim=c(0,150),breaks=seq(-0.15,0.1,by=0.01),xlab="Spearman's rho across samples")
abline(v=0.15)
hist(qualMat[,3],main=colnames(qualMat)[3],xlim=c(-0.15,0.1),ylim=c(0,150),breaks=seq(-0.15,0.1,by=0.01),xlab="Spearman's rho across samples")
abline(v=0.15)
hist(qualMat[,2],main=colnames(qualMat)[2],xlim=c(0,0.35),ylim=c(0,150),breaks=seq(0,0.35,by=0.01),xlab="Spearman's rho across samples")
abline(v=0.10)
hist(qualMat[,4],main=colnames(qualMat)[4],xlim=c(0,0.35),ylim=c(0,150),breaks=seq(0,0.35,by=0.01),xlab="Spearman's rho across samples")
abline(v=0.10)
dev.off()