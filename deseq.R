rm(list=ls())

my_file <- "projects/transcriptomics/homework/results2/Collibri_standard_protocol-HBR-Collibri-100_ng-2_S1_L001.featureCounts"

test <- read.csv(my_file, header=T, row.names="Geneid")

cts <- as.matrix(read.csv(my_file, header=TRUE, sep=""))


library(DESeq2)
library(ggplot2)




countFiles <- c("projects/transcriptomics/homework/results/Collibri_standard_protocol-HBR-Collibri-100_ng-2_S1_L001.featureCounts", 
                "projects/transcriptomics/homework/results/Collibri_standard_protocol-HBR-Collibri-100_ng-3_S2_L001.featureCounts", 
                "projects/transcriptomics/homework/results/Collibri_standard_protocol-UHRR-Collibri-100_ng-2_S3_L001.featureCounts",
                "projects/transcriptomics/homework/results/Collibri_standard_protocol-UHRR-Collibri-100_ng-3_S4_L001.featureCounts")

countDataList <- lapply(countFiles, function(file) {
  # Read the data from the file
  data <- read.table(file, header=TRUE, sep="", stringsAsFactors=FALSE)
  # Extract the count data from the last column
  counts <- data[,ncol(data)]
  # Return the count data as a named vector
  setNames(counts, data$Geneid)
})

# Combine the count data into a data frame
countData <- do.call(cbind, countDataList)

# Set the column names of the data frame
colnames(countData) <- c("HBR_ng2_s1", "HBR_ng3_s2", "UHRR_ng2_s3", "UHRR_ng3_s4")
colData <- data.frame(
  row.names = colnames(countData),
  condition = c("HBR", "HBR", "UHRR", "UHRR"),
  replicate = c(1, 2, 1, 2)
)

# Set the row names of the data frame to the gene IDs
rownames(countData) <- names(countDataList[[1]])
rownames(colData) <- colnames(countData)

condition <- factor(c("HBR", "HBR", "UHRR", "UHRR"))

condicolData <- data.frame(row.names=colnames(countData), condition)

countData <- countData[which(rowSums(countData) > 50),]

data <- as.data.frame(countData)

#data <- data[which(rowSums(data) > 50)]

# Create DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = colData,
                              design = ~condition)

dds <- DESeq(dds)
dds$condition <- factor(dds$condition, levels = c("HBR","UHRR"))
res <- results(dds, name = "condition_UHRR_vs_HBR")

resLFC <- lfcShrink(dds, coef="condition_UHRR_vs_HBR", type="apeglm")
resLFC

resOrdered <- res[order(res$pvalue),]

sum(res$padj < 0.05, na.rm=TRUE)
res05 <- results(dds, alpha=0.05)
summary(res05)
sum(res05$padj < 0.05, na.rm=TRUE)

library("IHW")
resIHW <- results(dds, filterFun=ihw)
summary(resIHW)
sum(resIHW$padj < 0.1, na.rm=TRUE)
metadata(resIHW)$ihwResult

plotMA(res, ylim=c(-2,2))
plotMA(resLFC, ylim=c(-2,2))

#idx <- identify(res$baseMean, res$log2FoldChange)
#rownames(res)[idx]

resNorm <- lfcShrink(dds, coef=2, type="normal")
resAsh <- lfcShrink(dds, coef=2, type="ashr")

par(mfrow=c(1,3), mar=c(4,4,2,1))
xlim <- c(1,1e5); ylim <- c(-3,3)
plotMA(resLFC, xlim=xlim, ylim=ylim, main="apeglm")
plotMA(resNorm, xlim=xlim, ylim=ylim, main="normal")
plotMA(resAsh, xlim=xlim, ylim=ylim, main="ashr")
par(mfrow=c(1,1))

d <- plotCounts(dds, gene=which.min(res$padj), intgroup="condition", 
                returnData=TRUE)
library("ggplot2")
ggplot(d, aes(x=condition, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400))
mcols(res)$description
write.csv(as.data.frame(resOrdered), 
          file="condition_treated_results.csv")
resSig <- subset(resOrdered, padj < 0.1)
resSig

vsd <- vst(dds, blind=FALSE)
vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)
head(assay(vsd), 3)
ntd <- normTransform(dds)
library("vsn")
meanSdPlot(assay(ntd))
meanSdPlot(assay(vsd))
meanSdPlot(assay(rld))

library("pheatmap")
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[,c("condition")])
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)
select <- order(rowMeans(counts(dds, normalized=TRUE)), decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[, c("condition")])
color_palette <- colorRampPalette(c("blue", "white", "red"))(50)
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)

pcaData <- plotPCA(vsd, intgroup=c("condition"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=condition, shape=condition)) +
  geom_point(size=2) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed(ratio = 1/3)

vsdata <-vst(dds, blind=FALSE)
vsdata1 <-varianceStabilizingTransformation(dds, blind=FALSE)

plotPCA(vsdata1)

plotDispEsts(dds)

res <- results(dds, contrast = c("condition", "HBR", "UHRR"))
res

plotMA(res, ylim=c(-2,2))

resLFC <- lfcShrink(dds, coef="condition HBR vs UHRR", type="apeglm")
resLFC

plotMA(res, ylim=c(-2,2))
plotMA(resLFC, ylim=c(-2,2))

resNorm <- lfcShrink(dds, coef=2, type="normal")
resAsh <- lfcShrink(dds, coef=2, type="ashr")

par(mfrow=c(1,3), mar=c(4,4,2,1))
xlim <- c(1,1e5); ylim <- c(-3,3)
plotMA(resLFC, xlim=xlim, ylim=ylim, main="apeglm")
plotMA(resNorm, xlim=xlim, ylim=ylim, main="normal")
plotMA(resAsh, xlim=xlim, ylim=ylim, main="ashr")

par(mfrow=c(1,1))
plotCounts(dds, gene=which.min(res$padj), intgroup="condition")

d <- plotCounts(dds, gene=which.min(res$padj), intgroup="condition", 
                returnData=TRUE)
library("ggplot2")
ggplot(d, aes(x=condition, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400))

mcols(res)$description

resOrdered <- res[order(res$pvalue),]
summary(res)

resSig <- subset(resOrdered, padj < 0.05)
resSig

vsd <- vst(dds, blind=FALSE)
vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)
head(assay(vsd), 3)

ntd <- normTransform(dds)
library("vsn")
meanSdPlot(assay(ntd))

meanSdPlot(assay(vsd))

meanSdPlot(assay(rld))

colData(dds)
ddsMF <- dds

levels(ddsMF$type)

levels(ddsMF$type) <- sub("-.*", "", levels(ddsMF$type))
levels(ddsMF$type)
design(ddsMF) <- formula(~ type + condition)
ddsMF <- DESeq(ddsMF)
resMF <- results(ddsMF)
head(resMF)
resMFType <- results(ddsMF,
                     contrast=c("type", "single", "paired"))
head(resMFType)


library("pheatmap")
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[,c("condition")])
rownames(df) <- colnames(ntd)
color_palette <- colorRampPalette(c("blue", "white", "red"))(50)
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)

pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)

pheatmap(assay(rld)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)

sampleDists <- dist(t(assay(vsd)))

library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$type, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

plotPCA(vsd, intgroup=c("condition"))

pcaData <- plotPCA(vsd, intgroup=c("condition"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=condition, shape=condition)) +
  geom_point(size=5) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed(ratio=8)

dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)
dds <- nbinomWaldTest(dds)

dds <- estimateSizeFactors(dds, controlGenes=ctrlGenes)
dds <- DESeq(dds)
# ------------------------------------------------------------------------------
# Get results
res <- results(dds)

file <- read.csv("projects/transcriptomics/homework/results/Collibri_standard_protocol-HBR-Collibri-100_ng-2_S1_L001.featureCounts", sep = "", header = TRUE, stringsAsFactors = FALSE)
test <- as.numeric(file[-1,7])
# Read count data
countDataList <- lapply(countFiles, function(x) {
  # Read count file
  print(x)
  df <- read.csv(x, sep = "", header = TRUE, stringsAsFactors = FALSE)
  # Extract count column
  counts <- as.numeric(df[-1,7])
  
  print(head(counts))
  counts
  #return(as.numeric(df[-1,6]))
})

# Combine count data into matrix
countData <- do.call(cbind, countDataList)

# Set column names
colnames(countData) <- c("HBR1", "HBR2", "UHRR1", "UHRR2")

# Create column data table
colData <- data.frame(condition = c("HBR", "HBR", "UHRR", "UHRR"))
rownames(colData) <- colnames(countData)

# Create DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = colData,
                              design = ~ condition)

# Run differential expression analysis
dds <- DESeq(dds)

# Get results
res <- results(dds)

# Filter results for adjusted p-value < 0.1
resSig <- subset(res, padj < 0.1)

# Order results by adjusted p-value
resSigOrdered <- resSig[order(resSig$padj),]

# Write results to file
write.csv(as.data.frame(resSigOrdered), file="DE_results.csv")

# Create volcano plot
plotMA(dds, ylim=c(-2,2))




sampleTable <- data.frame(
  condition = c("H", "U"),
  fileName = c("projects/transcriptomics/homework/results2/Collibri_standard_protocol-HBR-Collibri-100_ng-2_S1_L001.featureCounts", "projects/transcriptomics/homework/results2/Collibri_standard_protocol-UHRR-Collibri-100_ng-2_S3_L001.featureCounts")
)

countData <- lapply(sampleTable$fileName, read.table, header=TRUE, row.names=1)
countData <- do.call(cbind, countData)
colnames(countData) <- sampleTable$condition
dds <- DESeqDataSetFromMatrix(countData = countData, colData = sampleTable, design = ~condition)
dds <- DESeq(dds)
res <- results(dds)
resSig <- res[which(res$padj < 0.05),]
write.table(resSig, "DESeq2_results.txt", sep="\t", col.names=NA)
plot(log2(res$baseMean), -log10(res$padj), pch=20, cex=0.5, xlab="log2(counts)", ylab="-log10(padj)")
abline(h=-log10(0.05), col="red")
abline(v=c(-1, 1), col="blue")

