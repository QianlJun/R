###在满足Mean UMI counts > 1.0的基因中，按照log2FC从大到小排序选取top 20的基因作为该细胞亚群的候选特征基因
zhengchang_cluster_path<-"F:/zuoye/wanglab/data/Chromium-Singlecell-17-1409-01-人-4/zhengchang/result/4.cluster_specific_genes"
zhongliu_cluster_path<-"F:/zuoye/wanglab/data/Chromium-Singlecell-17-1409-01-人-4/zhongliu/result/4.cluster_specific_genes"
for (i in 2:10){
  cluster_path<-paste(zhongliu_cluster_path,"/kmeans_",i,sep="")
  setwd(cluster_path)
  expressionfile<-paste("kmeans_",i,"_clusters.differential_expression.csv",sep="")
  kmeans_expression<-read.csv(expressionfile,header = TRUE)
  for (j in 1:i ){
    kmeans_clusterj<-kmeans_expression[,c(2,3+(j-1)*3,4+(j-1)*3,5+(j-1)*3)]
    ##因子型变量先转换为字符型再转换为数值型
    kmeans_clusterj_UMI<-kmeans_clusterj[which(as.numeric(as.character(kmeans_clusterj[,2]))>1),]
    ##将log2FC从大到小排列
    kmeans_clusterj_UMI_fd<-kmeans_clusterj_UMI[order(kmeans_clusterj_UMI[,3],decreasing=T),]
    ##取top20的基因作为候选特征基因
    kmeans_clusterj_UMI_fd20<-kmeans_clusterj_UMI_fd[1:20,]
    kmeans_clusterj_UMI_fd20.gene<-kmeans_clusterj_UMI_fd20[,1]
    Newfile<-paste("Cluster_",j,"_symbol.csv",sep="")
    Newfile2<-paste("Cluster_",j,"-symbol.txt",sep="")
    write.csv(kmeans_clusterj_UMI_fd20,file=Newfile)
    write.table(kmeans_clusterj_UMI_fd20.gene,file=Newfile2,quote=FALSE,row.names=FALSE,col.names = FALSE)
  }
}

