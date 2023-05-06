rm(list=ls())


library(DESeq2)
library(ggplot2)



countFiles <- c("/home/jaro/projects/transcriptomics/homework/results/Collibri_standard_protocol-HBR-Collibri-100_ng-2_S1_L001.featureCounts", 
                "/home/jaro/projects/transcriptomics/homework/results/Collibri_standard_protocol-HBR-Collibri-100_ng-3_S2_L001.featureCounts", 
                "/home/jaro/projects/transcriptomics/homework/results/Collibri_standard_protocol-UHRR-Collibri-100_ng-2_S3_L001.featureCounts",
                "/home/jaro/projects/transcriptomics/homework/results/Collibri_standard_protocol-UHRR-Collibri-100_ng-3_S4_L001.featureCounts")

countDataList <- lapply(countFiles, function(file) {
  # Read the data from the file
  data <- read.table(file, header=TRUE, sep="", stringsAsFactors=FALSE)
  # Extract the count data from the last column
  counts <- data[,ncol(data)]
  # Return the count data as a named vector
  setNames(counts, data$Geneid)
})



countFilesKAPA <- c("/home/jaro/projects/transcriptomics/homework/results/KAPA_mRNA_HyperPrep_-HBR-KAPA-100_ng_total_RNA-2_S5_L001.featureCounts", 
                "/home/jaro/projects/transcriptomics/homework/results/KAPA_mRNA_HyperPrep_-HBR-KAPA-100_ng_total_RNA-3_S6_L001.featureCounts", 
                "/home/jaro/projects/transcriptomics/homework/results/KAPA_mRNA_HyperPrep_-UHRR-KAPA-100_ng_total_RNA-2_S7_L001.featureCounts",
                "/home/jaro/projects/transcriptomics/homework/results/KAPA_mRNA_HyperPrep_-UHRR-KAPA-100_ng_total_RNA-3_S8_L001.featureCounts")

countDataListKAPA <- lapply(countFilesKAPA, function(file) {
  # Read the data from the file
  data <- read.table(file, header=TRUE, sep="", stringsAsFactors=FALSE)
  # Extract the count data from the last column
  counts <- data[,ncol(data)]
  # Return the count data as a named vector
  setNames(counts, data$Geneid)
})



countData <- do.call(cbind, countDataList)

# Set the column names of the data frame
colnames(countData) <- c("HBR_ng2_s1", "HBR_ng3_s2", "UHRR_ng2_s3", "UHRR_ng3_s4")
colData <- data.frame(
  row.names = colnames(countData),
  condition = c("HBR", "HBR", "UHRR", "UHRR"),
  replicate = c(1, 2, 1, 2)
)



countDataKAPA <- do.call(cbind, countDataListKAPA)

# Set the column names of the data frame
colnames(countDataKAPA) <- c("HBR_s5", "HBR_s6", "UHRR_s7", "UHRR_s8")
colDataKAPA <- data.frame(
  row.names = colnames(countDataKAPA),
  condition = c("HBR", "HBR", "UHRR", "UHRR"),
  replicate = c(1, 2, 1, 2)
)




rownames(countData) <- names(countDataList[[1]])
rownames(colData) <- colnames(countData)

condition <- factor(c("HBR", "HBR", "UHRR", "UHRR"))

condicolData <- data.frame(row.names=colnames(countData), condition)

#countData <- countData[which(rowSums(countData) > 50),]

data <- as.data.frame(countData)



rownames(countDataKAPA) <- names(countDataListKAPA[[1]])
rownames(colDataKAPA) <- colnames(countDataKAPA)

condition <- factor(c("HBR", "HBR", "UHRR", "UHRR"))

condicolDataKAPA <- data.frame(row.names=colnames(countDataKAPA), condition)

#countDataKAPA <- countDataKAPA[which(rowSums(countDataKAPA) > 50),]

dataKAPA <- as.data.frame(countDataKAPA)



dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = colData,
                              design = ~condition)

dds <- DESeq(dds)
dds$condition <- factor(dds$condition, levels = c("HBR","UHRR"))
res <- results(dds, name = "condition_UHRR_vs_HBR")


ddsKAPA <- DESeqDataSetFromMatrix(countData = countDataKAPA,
                              colData = colDataKAPA,
                              design = ~condition)

ddsKAPA <- DESeq(ddsKAPA)
ddsKAPA$condition <- factor(ddsKAPA$condition, levels = c("HBR","UHRR"))
resKAPA <- results(ddsKAPA, name = "condition_UHRR_vs_HBR")





resLFC <- lfcShrink(dds, coef="condition_UHRR_vs_HBR", type="apeglm")
resLFC

resOrdered <- res[order(res$pvalue),]

sum(res$padj < 0.05, na.rm=TRUE)
res05 <- results(dds, alpha=0.05)
summary(res05)
sum(res05$padj < 0.05, na.rm=TRUE)




resLFCKAPA <- lfcShrink(ddsKAPA, coef="condition_UHRR_vs_HBR", type="apeglm")
resLFCKAPA

resOrderedKAPA <- resKAPA[order(resKAPA$pvalue),]

sum(resKAPA$padj < 0.05, na.rm=TRUE)
res05KAPA <- results(ddsKAPA, alpha=0.05)
summary(res05KAPA)
sum(res05KAPA$padj < 0.05, na.rm=TRUE)



library("IHW")


{r collibri-ma-plot}
resIHW <- results(dds, filterFun=ihw)
summary(resIHW)
sum(resIHW$padj < 0.1, na.rm=TRUE)
metadata(resIHW)$ihwResult

plotMA(res, ylim=c(-2,2))
plotMA(resLFC, ylim=c(-2,2))


{r kapa-ma-plot}
resIHWKAPA <- results(ddsKAPA, filterFun=ihw)
summary(resIHWKAPA)
sum(resIHWKAPA$padj < 0.1, na.rm=TRUE)
metadata(resIHWKAPA)$ihwResult

plotMA(resKAPA, ylim=c(-2,2))
plotMA(resLFCKAPA, ylim=c(-2,2))


{r collibri-ma-plot-3}
resNorm <- lfcShrink(dds, coef=2, type="normal")
resAsh <- lfcShrink(dds, coef=2, type="ashr")

par(mfrow=c(1,3), mar=c(4,4,2,1))
xlim <- c(1,1e5); ylim <- c(-3,3)
plotMA(resLFC, xlim=xlim, ylim=ylim, main="apeglm")
plotMA(resNorm, xlim=xlim, ylim=ylim, main="normal")
plotMA(resAsh, xlim=xlim, ylim=ylim, main="ashr")
par(mfrow=c(1,1))


{r kapa-ma-plot-3}
resNormKAPA <- lfcShrink(ddsKAPA, coef=2, type="normal")
resAshKAPA <- lfcShrink(ddsKAPA, coef=2, type="ashr")

par(mfrow=c(1,3), mar=c(4,4,2,1))
xlim <- c(1,1e5); ylim <- c(-3,3)
plotMA(resLFCKAPA, xlim=xlim, ylim=ylim, main="apeglm")
plotMA(resNormKAPA, xlim=xlim, ylim=ylim, main="normal")
plotMA(resAshKAPA, xlim=xlim, ylim=ylim, main="ashr")
par(mfrow=c(1,1))



library("ggplot2")



{r collibri-results}
d <- plotCounts(dds, gene=which.min(res$padj), intgroup="condition", 
                returnData=TRUE)
ggplot(d, aes(x=condition, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400))
mcols(res)$description
write.csv(as.data.frame(resOrdered), 
          file="condition_treated_results.csv")
resSig <- subset(resOrdered, padj < 0.1)
resSig


{r kapa-results}
dKAPA <- plotCounts(ddsKAPA, gene=which.min(resKAPA$padj), intgroup="condition", 
                returnData=TRUE)
ggplot(dKAPA, aes(x=condition, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400))
mcols(resKAPA)$description
write.csv(as.data.frame(resOrderedKAPA), 
          file="condition_treated_results_kapa.csv")
resSigKAPA <- subset(resOrderedKAPA, padj < 0.1)
resSigKAPA



{r collibri-mean-sd-plots}
vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)
head(assay(vsd), 3)
ntd <- normTransform(dds)
library("vsn")
meanSdPlot(assay(ntd))
meanSdPlot(assay(vsd))
meanSdPlot(assay(rld))


{r kapa-mean-sd-plots}
vsdKAPA <- varianceStabilizingTransformation(ddsKAPA, blind=FALSE)
rldKAPA <- rlog(ddsKAPA, blind=FALSE)
head(assay(vsdKAPA), 3)
ntdKAPA <- normTransform(ddsKAPA)

meanSdPlot(assay(ntdKAPA))
meanSdPlot(assay(vsdKAPA))
meanSdPlot(assay(rldKAPA))



library("pheatmap")


{r collibri-normal-heatmap}
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[,c("condition")])
rownames(df) <- colnames(ntd)
color_palette <- colorRampPalette(c("blue", "white", "red"))(50)
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)


{r collibri-vsd-heatmap}
pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)



{r collibri-rld-heatmap}
pheatmap(assay(rld)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)



sampleDists <- dist(t(assay(vsd)))


{r kapa-normal-heatmap}
selectKAPA <- order(rowMeans(counts(ddsKAPA,normalized=TRUE)),
                decreasing=TRUE)[1:20]
dfKAPA <- as.data.frame(colData(ddsKAPA)[,c("condition")])
rownames(dfKAPA) <- colnames(ntdKAPA)
color_palette <- colorRampPalette(c("blue", "white", "red"))(50)
pheatmap(assay(ntdKAPA)[selectKAPA,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=dfKAPA)


{r kapa-vsd-heatmap}
pheatmap(assay(vsdKAPA)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=dfKAPA)



{r kapa-rld-heatmap}
pheatmap(assay(rldKAPA)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=dfKAPA)



sampleDistsKAPA <- dist(t(assay(vsdKAPA)))


{r collibri-distmatrix}
library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$type, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)


{r kapa-distmatrix}
sampleDistMatrixKAPA <- as.matrix(sampleDistsKAPA)
rownames(sampleDistMatrixKAPA) <- paste(vsdKAPA$condition, vsdKAPA$type, sep="-")
colnames(sampleDistMatrixKAPA) <- NULL
colorsKAPA <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrixKAPA,
         clustering_distance_rows=sampleDistsKAPA,
         clustering_distance_cols=sampleDistsKAPA,
         col=colors)


{r collibri-pca}
pcaData <- plotPCA(vsd, intgroup=c("condition"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=condition, shape=condition)) +
  geom_point(size=5) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed(ratio=8)


{r kapa-pca}
pcaDataKAPA <- plotPCA(vsdKAPA, intgroup=c("condition"), returnData=TRUE)
percentVarKAPA <- round(100 * attr(pcaDataKAPA, "percentVar"))
ggplot(pcaDataKAPA, aes(PC1, PC2, color=condition, shape=condition)) +
  geom_point(size=5) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed(ratio=8)




library(EnhancedVolcano)


{r collibri-volkano}
  EnhancedVolcano(res,
    lab = rownames(res),
    x = 'log2FoldChange',
    y = 'pvalue',
    selectLab = c('VCAM1','KCTD12','ADAM12',
      'CXCL12','CACNB2','SPARCL1','DUSP1','SAMHD1','MAOA'),
    xlab = bquote(~Log[2]~ 'fold change'),
    pCutoff = 10e-14,
    FCcutoff = 2.0,
    pointSize = 4.0,
    labSize = 6.0,
    labCol = 'black',
    labFace = 'bold',
    boxedLabels = TRUE,
    colAlpha = 4/5,
    legendPosition = 'right',
    legendLabSize = 14,
    legendIconSize = 4.0,
    drawConnectors = TRUE,
    widthConnectors = 1.0,
    colConnectors = 'black')


{r kapa-volkano}
  EnhancedVolcano(resKAPA,
    lab = rownames(resKAPA),
    x = 'log2FoldChange',
    y = 'pvalue',
    xlab = bquote(~Log[2]~ 'fold change'),
    pCutoff = 10e-14,
    FCcutoff = 2.0,
    pointSize = 4.0,
    labSize = 2.0,
    labCol = 'black',
    labFace = 'bold',
    boxedLabels = TRUE,
    colAlpha = 4/5,
    legendPosition = 'right',
    legendLabSize = 14,
    legendIconSize = 4.0,
    drawConnectors = TRUE,
    widthConnectors = 1.0,
    colConnectors = 'black')



volcano_df <- data.frame(
  log2FC = res$log2FoldChange,
  Pvalue = res$pvalue,
  Gene = rownames(res)
)
EnhancedVolcano(
  volcano_df,
  lab = volcano_df$Gene,
  x = "log2FC",
  y = "Pvalue"
)




write.csv(as.data.frame(resLFC), file="/home/jaro/projects/transcriptomics/homework/results/res.csv", row.names = TRUE)

write.csv(as.data.frame(resLFCKAPA), file="/home/jaro/projects/transcriptomics/homework/results/resLFSKAPA.csv")



alpha <- 0.05 # Threshold on the adjusted p-value
cols <- densCols(res$log2FoldChange, -log10(res$pvalue))
plot.new()
plot(res$log2FoldChange, -log10(res$padj), col=cols, panel.first=grid(),
     main="Volcano plot", xlab="Effect size: log2(fold-change)", ylab="-log10(adjusted p-value)",
     pch=20, cex=0.6)

abline(v=0)
abline(v=c(-1,1), col="brown")
abline(h=-log10(alpha), col="brown")

gn.selected <- abs(res$log2FoldChange) > 2.5 & res$padj < alpha 
text(res$log2FoldChange[gn.selected],
     -log10(res$padj)[gn.selected],
     lab=rownames(res)[gn.selected ], cex=0.4)



