setwd("E:/D/wanglab/ATAC-seq/workdata/K")
setwd("C:/Users/MSI-MB/Desktop/workdata")
library("ChIPseeker")
library(clusterProfiler)
library("org.Hs.eg.db")
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
K<-readPeakFile("K_highquality_reproducible_fixed-width-peaks.bed")
peakAnno_K <- annotatePeak(K, tssRegion=c(-3000, 3000),
                         TxDb=txdb, annoDb="org.Hs.eg.db")
upsetplot(peakAnno_K)
W<-readPeakFile("W_highquality_reproducible_fixed-width-peaks.bed")
peakAnno_W <- annotatePeak(W, tssRegion=c(-3000, 3000),
                           TxDb=txdb, annoDb="org.Hs.eg.db")
upsetplot(peakAnno_W)
gene_K<-as.data.frame(peakAnno_K)$geneId
gene_W<-as.data.frame(peakAnno_W)$geneId
ego_K <- enrichGO(gene = gene_K, 
                OrgDb = org.Hs.eg.db, 
                ont = "CC", 
                pAdjustMethod = "BH", 
                qvalueCutoff = 0.05, 
                readable = TRUE)
dotplot(ego_K, showCategory=50)
ego_W <- enrichGO(gene = gene_W, 
                  OrgDb = org.Hs.eg.db, 
                  ont = "CC", 
                  pAdjustMethod = "BH", 
                  qvalueCutoff = 0.05, 
                  readable = TRUE)
dotplot(ego_W, showCategory=50)
peakAnno_K_data<-as.data.frame(peakAnno_K)
write.table(peakAnno_K_data,file="peakAnno_K_data.txt",sep = "\t",quote = FALSE,col.names = FALSE)
