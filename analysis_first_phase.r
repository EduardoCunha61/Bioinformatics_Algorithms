if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")

ianstall.packages('DESeq2')
install.packages("pheatmap")


library(DESeq2)
library(pheatmap)

#1- load the mapping data and create a DESeqDataSet object
# load the read count matrix
total<- read.table("C2C11C17.tab", h=T, row.names=1)
View(total)


# define the conditions to compare
condition <- factor(c("C2", "C2", "C2","C2", "C2", "C2", "C11", "C11", "C11", "C11", "C11", "C11", "C17", "C17", "C17","C17","C17","C17"))


# pre-filter low count genes before running the DESeq2 functions
total <- total[ rowSums(total) > 1, ]
View(total)
dim(total)

# define the conditions of the samples
cd2=data.frame(c("C2", "C2", "C2","C2", "C2", "C2", "C11", "C11", "C11", "C11", "C11", "C11", "C17", "C17", "C17","C17","C17","C17"))
colnames(cd2)[1]="condition"
rownames(cd2)=colnames(total)
View(cd2)

#2- Run Differential expression test
# create an object for Diff Expr analysis
dds2 <- DESeqDataSetFromMatrix(countData = total, colData = cd2, design = ~ condition)

# differential expression of samples that belong to C2 vs C11 vs C17
dds2 <- DESeq(dds2)
res2 <- results(dds2)
res2

# order our results table by the smallest adjusted p value:
resOrdered <- res2[order(res2$padj),]

# summarize some basic tallies using the summary function
summary(res2)

# many adjusted p-values were less than 0.1?
sum(res2$padj < 0.1, na.rm=TRUE)


#3- Exploring and exporting results
# MA-plot
plotMA(res2, main="DESeq2", ylim=c(-2,2))


# plot the read counts for the most significant gene
plotCounts(dds2, gene=which.min(res2$padj), intgroup="condition")


# Exporting results to CSV files
write.csv(as.data.frame(resOrdered),file="condition_cancer_results.csv")

#4- Data transformations and visualization
# VST: varianceStabilizingTransformation
vsd2 <- varianceStabilizingTransformation(dds2, blind=FALSE)
head(assay(vsd2), 3)

# compare previous command with head(counts(dds), 3)

# Heatmap of the count matrix
select <- rownames(head(resOrdered,20))
vsd.counts <- assay(vsd2)[select,]
df2 <- as.data.frame(colData(dds2)[,c("condition")])

# maintain row order
pheatmap(vsd.counts, cluster_rows=FALSE)

# cluster by row and column
pheatmap(vsd.counts)
