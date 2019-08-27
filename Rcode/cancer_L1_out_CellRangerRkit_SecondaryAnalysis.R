library(cellrangerRkit)
##Loading gene expression data for secondary analysis
cellranger_pipestance_path<-"F:/workdata/cancer_L1"
gbm <- load_cellranger_matrix(cellranger_pipestance_path)
analysis_results <- load_cellranger_analysis_results(cellranger_pipestance_path)
exprs(gbm) # expression matrix
fData(gbm) # data frame of genes
pData(gbm) # data frame of cell barcodes
##access the t-SNE projection and plot the cells colored by UMI counts
tsne_proj <- analysis_results$tsne
visualize_umi_counts(gbm,tsne_proj[c("TSNE.1","TSNE.2")],limits=c(3,5),marker_size=1.2)
##filter unexpressed genes,normalize the UMI counts for each barcode
##and use the log-transformed gene-barcode matrix
use_genes <- get_nonzero_genes(gbm)
gbm_bcnorm <- normalize_barcode_sums_to_median(gbm[use_genes,])
gbm_log <- log_gene_bc_matrix(gbm_bcnorm,base=10)
print(dim(gbm_log))
##Unbiased analysis using clustering results
##k-means clustering , k form 2 to 10
n_clu <- 2:10
km_res <- analysis_results$clustering # load pre-computed kmeans results
clu_res <- sapply(n_clu, function(x) km_res[[paste("kmeans",x,"clusters",sep="_")]]$Cluster)
colnames(clu_res) <- sapply(n_clu, function(x) paste("kmeans",x,sep="."))
visualize_clusters(clu_res,tsne_proj[c("TSNE.1","TSNE.2")])
##Analyzing cluster specific genes 
##identify gene markers that are specific to a particular cell population
example_K <- 5 # number of clusters (use "Set3" for brewer.pal below if example_K > 8)
example_col <- rev(brewer.pal(example_K,"Set2")) # customize plotting colors
cluster_result <- analysis_results$clustering[[paste("kmeans", example_K,"clusters",sep="_")]]
visualize_clusters(cluster_result$Cluster,tsne_proj[c("TSNE.1","TSNE.2")],colour=example_col)
##compare the mean expression between a class of cells and the remaining ones
##and then prioritize genes by how highly expressed they are in the class of interest
# sort the cells by the cluster labels
cells_to_plot <- order_cell_by_clusters(gbm, cluster_result$Cluster)
# order the genes from most up-regulated to most down-regulated in each cluster
prioritized_genes <- prioritize_top_genes(gbm, cluster_result$Cluster, "sseq", min_mean=0.5)
## output the top genes specific to each cluster to a local folder
# in this case,we output all the top 10 gene symbels for thee 5 cluster to file
output_folder <-"F:/workdata/cancer_L1/gene_sets"
write_cluster_specific_genes(prioritized_genes, output_folder, n_genes=10)
# create values and axis annotations for pheatmap
gbm_pheatmap(log_gene_bc_matrix(gbm), prioritized_genes, cells_to_plot,
             n_genes=4, colour=example_col, limits=c(-1,2))
##find out the size of each cluster
cell_composition(cluster_result$Cluster,
                 anno=c("cluster1","cluster2","cluster3","cluster4","cluster5"))
##Visualizing signatures of known gene markers
##simultaneously plot the expression of 6 gene markers
genes <- c("GABRP","CXCL17","LAPTM5","IGFBP7","DCN")
tsne_proj <- analysis_results$tsne
visualize_gene_markers(gbm_log,genes,tsne_proj[c("TSNE.1","TSNE.2")],limits=c(0,1.5))
