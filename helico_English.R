setwd("C:/Users/Dell/Documents/work/data/helico")
library("DESeq2")
library("RColorBrewer")
library("gplots")
library("ggplot2")
library("AnnotationDbi")
library("org.Hs.eg.db")
library(AnnotationForge)


#metadata
metadata <- read.csv("helico_meta.csv", row.names = "Sample")
metadata$Replicate <- as.factor(metadata$Replicate)
metadata$Cell_line <- as.factor(metadata$Cell_line)
metadata$Strain <- as.factor(metadata$Strain)
metadata$cell_strain <- paste(metadata$Cell_line, metadata$Strain)
metadata$cell_strain <- as.factor(metadata$cell_strain)


#counts
helico_counts <- as.matrix(read.csv("counts_R.csv",row.names="Gene"))

#removing 5 lines from the end
tail(helico_counts)
n <- dim(helico_counts)[1]
helico_counts <- helico_counts[1:(n-5),]
tail(helico_counts)

#checking compatibility
row.names(metadata)==colnames(helico_counts)

#removing genes without expression
keep <- rowSums(helico_counts) > 0
helico_counts <- helico_counts[keep,]
dim(helico_counts)

# Compiling Deseq Dataset without design
CellTable <- DESeqDataSetFromMatrix(
  countData = helico_counts,
  colData = metadata,
  design = ~ 1)

#Prefiltering
CellTable <-  estimateSizeFactors(CellTable)

#Selection of genes with abundat expressiion (2 normalized reads in at least 3 samples)
keep <- rowSums(counts(CellTable, normalized=T) >= 2) >= 3
CellTable <- CellTable[keep,]
dim(counts(CellTable))


#VST transformation for visualisation
vsd <- vst(CellTable)

#PCA plot and save in PDF 
pdf("PCA helico rakuliin.pdf")
plotPCA(vsd, intgroup=c("Strain"))
dev.off()

#Heatmap on 500 most variable genes
topVarGenes <- head( order( rowVars( assay(vsd) ), decreasing=TRUE ), 500 )

pdf("proov heatmap.pdf")
heatmap.2( assay(vsd)[ topVarGenes, ], scale="row",
           trace="none", dendrogram="column", labCol = c(as.character(vsd$cell_strain)), margins = c(6,5),
           col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255))
dev.off()

#-------IHH cell-line analysis-------

#Extracting the relevant cell-line
IHH <- CellTable[,CellTable$Cell_line==c("IHH")]
  
# Re-filtering 
IHH <-  estimateSizeFactors(IHH)

#Selecting high abundance genes
keep <- rowSums(counts(IHH, normalized=T) >= 2) >= 3
IHH <- IHH[keep,]
dim(counts(IHH)) 


#VST transformation for visualisation
vsd_IHH <- vst(IHH)

#PCA plot and saving PDF 
#For visualizing other parameters, change "intgroup"
pdf("PCA helico tyvi IHH.pdf")
plotPCA(vsd_IHH, intgroup=c("Strain"))
dev.off()

pdf("PCA helico replikaat IHH.pdf")
plotPCA(vsd_IHH, intgroup=c("Replicate"))
dev.off()

#Remove replicate information as batch effect
#PCA plot without batch effect. 
mat <- assay(vsd_IHH)
mat <- limma::removeBatchEffect(mat, vsd_IHH$Replicate)
assay(vsd_IHH) <- mat

pdf("PCA IHH bakterityvi replikaadi efekt eemaldatud.pdf")
plotPCA(vsd_IHH, intgroup=c("Strain"))
dev.off()

#DE analysis
design(IHH) <- ~Replicate+Strain

#Specify control group
IHH$Strain <- relevel(IHH$Strain, ref = "control")

#DEseq
IHH <- DESeq(IHH)
res <- results(IHH)
res

#See all groups
levels(IHH$Strain)


#Change comparison groups
res1 <-  results(IHH, contrast=c("Strain","J99","control"))
res2 <-  results(IHH, contrast=c("Strain","7.13_wt","control"))
res3 <-  results(IHH, contrast=c("Strain","X47", "7.13_wt"))
res4 <-  results(IHH, contrast=c("Strain","X47", "J99"))
res5 <-  results(IHH, contrast=c("Strain","J99", "7.13_wt"))

#Display the number of significant genes
sum( res5$padj < 0.05, na.rm=TRUE )

#order by FDR
resOrdered5 <- res5[order(res5$padj),]

#Save only significant results in a file
resSig5 <- subset(resOrdered5, padj < 0.05)
resSig5


#Annotation of tables
columns(org.Hs.eg.db)

#ens.str <- substr(rownames(res), 1, 15)
ens.str <- rownames(resSig5)
resSig5$symbol <- mapIds(org.Hs.eg.db,
                             keys=ens.str,
                             column="SYMBOL",
                             keytype="ENSEMBL",
                             multiVals="first")
resSig5$genename <- mapIds(org.Hs.eg.db,
                               keys=ens.str,
                               column="GENENAME",
                               keytype="ENSEMBL",
                               multiVals="first")


#Check that all required fields are present
resSig5

#Save the results
write.csv(resSig5, file="Significant padj genes J99 vs 7.13_wt Sept28_2021.csv")


#Plot the expression of single genes
#plot counts
pdf("TRAF1 in IHH.pdf", width=10)
plotCounts(IHH, gene="ENSG00000056558", intgroup="Strain", normalized=T, xlab= "Strain",main="TRAF1 in IHH")
dev.off()


pdf("CD44 in IHH.pdf", width=10)
plotCounts(IHH, gene="ENSG00000026508", intgroup="Strain", normalized=T, xlab= "Strain",main="CD44 in IHH")
dev.off()

