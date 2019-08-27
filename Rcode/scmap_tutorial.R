setwd("E:/D/wanglab/Single-cell/Linkai/scmap")
library(SingleCellExperiment)
library(scmap)
head(ann,10)
write.table(ann,file="ann.txt",sep="\t",quote = F,col.names = T,row.names = T)
write.table(yan,file="yan.txt",sep = "\t",quote = F,row.names = T,col.names = T)
yan[1:3,1:3]

sce <- SingleCellExperiment(assays = list(normcounts = as.matrix(yan)), colData = ann)
logcounts(sce) <- log2(normcounts(sce) + 1)
# use gene names as feature symbols
rowData(sce)$feature_symbol <- rownames(sce)
isSpike(sce, "ERCC") <- grepl("^ERCC-", rownames(sce))
# remove features with duplicated names
sce <- sce[!duplicated(rownames(sce)), ]
sce

## Feature selection
sce <- selectFeatures(sce, suppress_plot = FALSE,n_features = 500)
table(rowData(sce)$scmap_features)
sce <- indexCluster(sce)
head(metadata(sce)$scmap_cluster_index)
heatmap(as.matrix(metadata(sce)$scmap_cluster_index))
