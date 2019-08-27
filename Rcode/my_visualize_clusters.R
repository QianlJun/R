my.visualize_clusters<-function (cluster_result, projection, colour = NULL, alpha = 1, 
                                 marker_size = 0.1, title = NULL, legend_anno = NULL) 
{
  cluster_ids <- sort(unique(as.vector(cluster_result)))
  if (is.vector(cluster_result)) {
    if (dim(projection)[1] != length(cluster_result)) 
      stop("The number of labels and the number of projected points do not match!\n")
  }
  else {
    if (dim(projection)[1] != dim(cluster_result)[1]) {
      stop("The number of labels and the number of projected points do not match!\n")
    }
  }
  if (!is.null(legend_anno)) {
    if (length(legend_anno) != length(cluster_ids)) 
      stop("The length of legend_anno is not the same as the number of unique labels!\n")
  }
  if (!is.null(colour)) {
    if (length(colour) != length(cluster_ids)) 
      stop("The length of colour is not the same as the number of unique labels!\n")
    names(colour) <- cluster_ids
    if (!is.null(legend_anno)) {
      names(colour) <- legend_anno
    }
  }
  projection_names <- colnames(projection)
  colnames(projection) <- c("Component.1", "Component.2")
  proj_clu <- data.frame(cbind(projection, cluster_result))
  proj_clu_melt <- melt(proj_clu, id.vars = c("Component.1", 
                                              "Component.2"))
  if (!is.null(legend_anno)) {
    cid_idx <- unlist(lapply(cluster_result, function(x) which(cluster_ids %in% 
                                                                 x)))
    proj_clu_melt$value <- factor(legend_anno[cid_idx])
  }
  else {
    proj_clu_melt$value <- factor(proj_clu_melt$value)
  }
  p <- ggplot(proj_clu_melt, aes(Component.1, Component.2)) + 
    geom_point(aes(colour = value), size = marker_size, 
               alpha = alpha) + facet_wrap(~variable) + guides(col = guide_legend(title = "Hello", 
                                                                                  override.aes = list(size = 3))) + labs(x = projection_names[1], 
                                                                                                                         y = projection_names[2]) + theme_bw()
  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }
  if ((!is.null(colour))) {
    p <- p + scale_color_manual(values = colour)
  }
  if (is.vector(cluster_result)) {
    p <- p + theme(plot.title = element_text(hjust = 0.5), 
                   legend.key = element_blank(), strip.text.x = element_blank(), 
                   panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  }
  else {
    p <- p + theme(plot.title = element_text(hjust = 0.5), 
                   legend.key = element_blank(), panel.grid.major = element_blank(), 
                   panel.grid.minor = element_blank())
  }
  p<-p+theme(strip.background=element_rect(fill="white", colour = "black", size=rel(1), linetype = 1),strip.text=element_text(size=rel(2.0)))
  return(p)
}
library(reshape2)
my.visualize_clusters(clu_res,tsne_proj[c("TSNE.1","TSNE.2")])
