## 1.Setup the Seurat Object
library(dplyr)
library(Seurat)
# Load the PBMC dataset
pbmc.data<-Read10X(data.dir="E:/D/wanglab/Single-cell/Linkai/filtered_gene_bc_matrices/hg19")
# Examine the memory savings between regular and sparse matrices
dense.size<-object.size(x=as.matrix(x=pbmc.data))
dense.size
sparse.size <- object.size(x = pbmc.data)
sparse.size
dense.size / sparse.size
# Initialize the Seurat object with the raw (non-normalized data).
# Keep all features expressed in >= 3 cells (~0.1% of the data). Keep all cells with at least 200 detected features
pbmc <- CreateSeuratObject(counts = pbmc.data, min.cells = 3, min.features = 200, project = "10X_PBMC")
pbmc
## 2.Standard pre-processing workflow
## QC and selecting cells for further analysis
# The number of features and UMIs (nFeature_RNA and nCount_RNA) are automatically calculated for every object by Seurat.
# For non-UMI data, nCount_RNA represents the sum of the non-normalized values within a cell
# We calculate the percentage of mitochondrial features here and store it in object metadata as `percent.mito`.
# We use raw count data since this represents non-transformed and non-log-normalized counts
# The % of UMI mapping to MT-features is a common scRNA-seq QC metric.
mito.features <- grep(pattern = "^MT-", x = rownames(x = pbmc), value = TRUE)
percent.mito <- Matrix::colSums(x = GetAssayData(object = pbmc, slot = 'counts')[mito.features, ]) / Matrix::colSums(x = GetAssayData(object = pbmc, slot = 'counts'))

# The [[ operator can add columns to object metadata, and is a great place to stash QC stats
pbmc[['percent.mito']] <- percent.mito
VlnPlot(object = pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used for anything 
# calculated by the object, i.e. columns in object metadata, PC scores etc.
# Since there is a rare subset of cells with an outlier level of high mitochondrial percentage
# and also low UMI content, we filter these as well
FeatureScatter(object = pbmc, feature1 = "nCount_RNA", feature2 = "percent.mito")
FeatureScatter(object = pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

# We filter out cells that have unique feature counts over 2,500 or less than 200
pbmc <- subset(x = pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mito < 0.05)

## 3.Normalizing the data
pbmc <- NormalizeData(object = pbmc,normalization.method = "LogNormalize", scale.factor = 1e4)

## 4.Detection of variable features across the single cells
pbmc <- FindVariableFeatures(object = pbmc, selection.method = 'mean.var.plot', mean.cutoff = c(0.0125, 3), dispersion.cutoff = c(0.5, Inf))
length(x = VariableFeatures(object = pbmc))

## 5.Scaling the data and removing unwanted sources of variation
pbmc <- ScaleData(object = pbmc, features = rownames(x = pbmc), vars.to.regress = c("nCount_RNA", "percent.mito"))
features = rownames(x = pbmc)
length(features)
str(pbmc)
                    
## 6.Perform linear dimensional reduction
pbmc <- RunPCA(object = pbmc, features = VariableFeatures(object = pbmc), verbose = FALSE)
# Examine and visualize PCA results a few different ways
print(x = pbmc[['pca']], dims = 1:5, nfeatures = 5, projected = FALSE)
VizDimLoadings(object = pbmc, dims = 1:2)
DimPlot(object = pbmc)
# ProjectDim scores each feature in the dataset (including features not included in the PCA) based on their correlation 
# with the calculated components. Though we don't use this further here, it can be used to identify markers that 
# are strongly correlated with cellular heterogeneity, but may not have passed through variable feature selection. 
# The results of the projected PCA can be explored by setting `projected = TRUE`in the functions above
pbmc <- ProjectDim(object = pbmc)
DimHeatmap(object = pbmc, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(object = pbmc, dims = 1:12, cells = 500, balanced = TRUE)

## 7.Determine statistically significant principal components
# NOTE: This process can take a long time for big datasets, comment out for expediency.
# More approximate techniques such as those implemented in ElbowPlot() can be used to reduce computation time
pbmc <- JackStraw(object = pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(object = pbmc, dims = 1:20)
JackStrawPlot(object = pbmc, dims = 1:12)
ElbowPlot(object = pbmc)

## 8.Cluster the cells
pbmc <- FindNeighbors(object = pbmc, dims = 1:10)
pbmc <- FindClusters(object = pbmc, resolution = 0.6)

## 9.Run Non-linear dimensional reduction (tSNE)
pbmc <- RunTSNE(object = pbmc, dims = 1:10)
# note that you can set `label = TRUE` or use the LabelClusters function to help label individual clusters
DimPlot(object = pbmc, reduction = 'tsne')
saveRDS(pbmc, file = "E:/D/wanglab/Single-cell/Linkai/filtered_gene_bc_matrices/hg19/pbmc_tutorial.rds")

## 10.Finding differentially expressed features (cluster biomarkers)
# find all markers of cluster 1
cluster1.markers <- FindMarkers(object = pbmc, ident.1 = 1, min.pct = 0.25)
head(x = cluster1.markers, n = 5)
# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(object = pbmc, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
head(x = cluster5.markers, n = 5)
# find markers for every cluster compared to all remaining cells, report only the positive ones
pbmc.markers <- FindAllMarkers(object = pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
cluster1.markers <- FindMarkers(object = pbmc, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
VlnPlot(object = pbmc, features = c("MS4A1", "CD79A"))
# you can plot raw counts as well
VlnPlot(object = pbmc, features = c("NKG7", "PF4"), slot = 'counts', log = TRUE)
FeaturePlot(object = pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", "CD8A"))
pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC) -> top10
DoHeatmap(object = pbmc, features = top10$gene) + NoLegend()

## 11.Assigning cell type identity to clusters
new.cluster.ids <- c("CD4 T cells", "CD14+ Monocytes", "B cells", "CD8 T cells", "FCGR3A+ Monocytes", "NK cells", "Dendritic cells", "Megakaryocytes")
names(x = new.cluster.ids) <- levels(x = pbmc)
pbmc <- RenameIdents(object = pbmc, new.cluster.ids)
DimPlot(object = pbmc, reduction = 'tsne', label = TRUE, pt.size = 0.5) + NoLegend()

## 12.Further subdivisions within cell types
# First lets stash our identities for later
pbmc <- StashIdent(object = pbmc, save.name = "ClusterNames_0.6")
# Note that you don't need to recalculate the SNN, and can simply put: 
pbmc <- FindClusters(object = pbmc, dims = 1:10, resolution = 0.8)
# Demonstration of how to plot two tSNE plots side by side, and how to color points based on different criteria
plot1 <- DimPlot(object = pbmc, label = TRUE)
plot2 <- DimPlot(object = pbmc, group.by = "ClusterNames_0.6", label = TRUE)
CombinePlots(plots = list(plot1, plot2), legend = 'none')
# Find discriminating markers
tcell.markers <- FindMarkers(object = pbmc, ident.1 = 0, ident.2 = 1)
head(x = tcell.markers)
# Most of the markers tend to be expressed in C1 (i.e. S100A4). However, we can see that CCR7 is upregulated in 
# C0, strongly indicating that we can differentiate memory from naive CD4 cells.
# `cols` demarcates the color palette from low to high expression
FeaturePlot(object = pbmc, features = c("S100A4", "CCR7"), cols = c("green", "blue"))
Idents(object = pbmc) <- 'ClusterNames_0.6'
saveRDS(pbmc, file = "E:/D/wanglab/Single-cell/Linkai/filtered_gene_bc_matrices/hg19/pbmc3k_final.rds")
