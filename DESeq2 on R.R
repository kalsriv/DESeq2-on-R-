library(DESeq2)
library(edgeR)
library(limma)
library(Biobase)
library(dplyr)
library(ggplot2)
library(org.Hs.eg.db)
library(pheatmap)

load("/home/srivastava/Documents/HTSeq_Example/gilad_eset.RData")
data <- gilad.eset
head(data)
print(data)
mtData <-pData(data)
cpm.mat <- log(cpm(exprs(data)))
mean.vec <- apply(cpm.mat, 1, mean)
sdvec <- apply(cpm.mat, 1, sd)
plot(mean.vec, sdvec, pch =".", ylab="sd", xlab="Average logCPM")

colnames(data)
dds <- DESeqDataSetFromMatrix(countData = exprs(data), colData = mtData, design = ~gender)
dds <- DESeq(dds)
plotDispEsts(dds)


# need transformed values to generate a heat plot of high differential expression
vsd <- vst(dds, blind=FALSE) # Variance stabilizing transformation
rld <- rlog(dds, blind=FALSE) # Controls for rows with few counts
ntd <- normTransform(dds) # Transforms data for plotting
top20 <- order(rowMeans(counts(dds,normalized=TRUE)),
               decreasing=TRUE)[1:20] # Average counts for each gene and Select top 20 most abundant genes
pheatmap(assay(ntd)[top20,])




#Next we would work with the gene names and try to find significant and positive log fold genes
ddsDF <- results(dds)
class(ddsDF) #This is not a data frame yet
ddsDF <- as.data.frame(ddsDF) #This is a data frame
head(ddsDF)
gene.list <- rownames(ddsDF)
head(gene.list)

#Let us add the gene names from the annotations

geneinfo <- AnnotationDbi::select(org.Hs.eg.db, keys=rownames(dds),
                                  columns=c("ENSEMBL","SYMBOL","GENENAME"), 
                                  keytype="ENSEMBL")
geneinfo %>% head()#... since this is huge number of entries

# Let us work only with significant values and generate the list of genes that are significant
ddsDF05 <- ddsDF[which(ddsDF$pvalue <= 0.05),]
geneinfo05 <- AnnotationDbi::select(org.Hs.eg.db, keys=rownames(ddsDF05),
                                  columns=c("ENSEMBL","SYMBOL","GENENAME"), 
                                  keytype="ENSEMBL")

ddsDF05_positive <- ddsDF05[ddsDF05$log2FoldChange >= 0,]
geneinfo05_positive <- AnnotationDbi::select(org.Hs.eg.db, keys=rownames(ddsDF05_positive),
                                    columns=c("ENSEMBL","SYMBOL","GENENAME"), 
                                    keytype="ENSEMBL")

geneinfo05_positive
#This would generate a list of genes that are significant and have log2fold change as positive


#Now some PCA analysis

exprsData <-  as.data.frame(exprs(data))
exprsMat <- as.matrix(exprsData)
df_pca <- prcomp(exprsMat)
plot(df_pca$x[,1], df_pca$x[,2])
#plot(df_pca$x[,3], df_pca$x[,4])

