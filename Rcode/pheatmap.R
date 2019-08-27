install.packages("pheatmap")
library(pheatmap)
setwd("E:/D/wanglab/hMeDIP-seq/workdata")
# 构建测试数据集
test = matrix(rnorm(200), 20, 10)
test[1:10, seq(1, 10, 2)] = test[1:10, seq(1, 10, 2)] + 3
test[11:20, seq(2, 10, 2)] = test[11:20, seq(2, 10, 2)] + 2
test[15:20, seq(2, 10, 2)] = test[15:20, seq(2, 10, 2)] + 4
colnames(test) = paste("Test", 1:10, sep = "")
rownames(test) = paste("Gene", 1:20, sep = "")
head(test)
# 默认绘图
pheatmap(test)
# scale = "row"参数对行进行归一化
pheatmap(test, scale = "row")
# clustering_method参数设定不同聚类方法，默认为"complete",可以设定为'ward', 'ward.D', 'ward.D2', 'single', 'complete', 'average', 'mcquitty', 'median' or 'centroid'
pheatmap(test,scale = "row", clustering_method = "average")
# clustering_distance_rows = "correlation"参数设定行聚类距离方法为Pearson corralation，默认为欧氏距离"euclidean"
pheatmap(test, scale = "row", clustering_distance_rows = "correlation")
# color参数自定义颜色
pheatmap(test, color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
# cluster_row = FALSE参数设定不对行进行聚类
pheatmap(test, cluster_row = FALSE)
# legend_breaks参数设定图例显示范围，legend_labels参数添加图例标签
pheatmap(test, legend_breaks = c(1:5), legend_labels = c("1.0","2.0","3.0","4.0","5.0"))
# legend = FALSE参数去掉图例
pheatmap(test, legend = FALSE)
# border_color参数设定每个热图格子的边框色
pheatmap(test, border_color = "red")
# border=FALSE参数去掉边框线
pheatmap(test, border=FALSE)
# show_rownames和show_colnames参数设定是否显示行名和列名
pheatmap(test,show_rownames=F,show_colnames=F)
# treeheight_row和treeheight_col参数设定行和列聚类树的高度，默认为50
pheatmap(test, treeheight_row = 30, treeheight_col = 50)
# display_numbers = TRUE参数设定在每个热图格子中显示相应的数值，number_color参数设置数值字体的颜色
pheatmap(test, display_numbers = TRUE,number_color = "blue")
# number_format = "%.1e"参数设定数值的显示格式
pheatmap(test, display_numbers = TRUE, number_format = "%.1e")
# 自定义数值的显示方式
pheatmap(test, display_numbers = matrix(ifelse(test > 5, "*", ""), nrow(test)))
# cellwidth和cellheight参数设定每个热图格子的宽度和高度，main参数添加主标题
pheatmap(test, cellwidth = 15, cellheight = 12, main = "Example heatmap")
image.png

# 构建列注释信息
annotation_col = data.frame(
  CellType = factor(rep(c("CT1", "CT2"), 5)), 
  Time = 1:5
)
rownames(annotation_col) = paste("Test", 1:10, sep = "")
head(annotation_col)
# 构建行注释信息
annotation_row = data.frame(
  GeneClass = factor(rep(c("Path1", "Path2", "Path3"), c(10, 4, 6)))
)
rownames(annotation_row) = paste("Gene", 1:20, sep = "")
head(annotation_row)
# annotation_col参数添加列注释信息
pheatmap(test, annotation_col = annotation_col)
# annotation_legend = FALSE参数去掉注释图例
pheatmap(test, annotation_col = annotation_col, annotation_legend = FALSE)
# annotation_col和annotation_row参数同时添加行和列的注释信息
pheatmap(test, annotation_row = annotation_row, annotation_col = annotation_col)
# 自定注释信息的颜色列表
ann_colors = list(
  Time = c("white", "firebrick"),
  CellType = c(CT1 = "#1B9E77", CT2 = "#D95F02"),
  GeneClass = c(Path1 = "#7570B3", Path2 = "#E7298A", Path3 = "#66A61E")
)
head(ann_colors)
# annotation_colors设定注释信息的颜色
pheatmap(test, annotation_col = annotation_col, annotation_colors = ann_colors, main = "Title")
pheatmap(test, annotation_col = annotation_col, annotation_row = annotation_row, 
         annotation_colors = ann_colors)
pheatmap(test, annotation_col = annotation_col, annotation_colors = ann_colors[2])
image.png

# gaps_row = c(10, 14)参数在第10和14行处添加gap, 要求对行不进行聚类
pheatmap(test, annotation_col = annotation_col, cluster_rows = FALSE, gaps_row = c(10, 14))
# cutree_col = 2参数将列按聚类树的结果分成两部分, 要求对列进行聚类
pheatmap(test, annotation_col = annotation_col, cluster_rows = FALSE, gaps_row = c(10, 14),
         cutree_col = 2)
# 对行和列都不聚类，自定义划分行和列的gap
pheatmap(test, annotation_col = annotation_col, cluster_rows = FALSE, cluster_cols = FALSE, 
         gaps_row = c(6, 10, 14), gaps_col = c(2, 5, 8))
# 自定义行的标签名
labels_row = c("", "", "", "", "", "", "", "", "", "", "", "", "", "", "", 
               "", "", "Il10", "Il15", "Il1b")
# labels_row参数添加行标签
pheatmap(test, annotation_col = annotation_col, labels_row = labels_row)
# 自定义聚类的距离方法
drows = dist(test, method = "minkowski")
dcols = dist(t(test), method = "minkowski")
# clustering_distance_rows和clustering_distance_cols参数设定行和列的聚类距离方法
pheatmap(test, clustering_distance_rows = drows, clustering_distance_cols = dcols)
# fontsize参数设定标签字体大小，filename参数设定图片保存名称
pheatmap(test, cellwidth = 15, cellheight = 12, fontsize = 8, filename = "test.pdf")
aa=pheatmap(test,scale="row")  #热图，归一化，并聚类
# 简要查看热图对象的信息
summary(aa)
order_row = aa$tree_row$order  #记录热图的行排序
order_col = aa$tree_col$order    #记录热图的列排序
datat = data.frame(test[order_row,order_col])   # 按照热图的顺序，重新排原始数据
datat = data.frame(rownames(datat),datat,check.names =F)  # 将行名加到表格数据中
colnames(datat)[1] = "geneid" 
write.table(datat,file="reorder.txt",row.names=FALSE,quote = FALSE,sep='\t')  #输出结果，按照热图中的顺序
sessionInfo()
