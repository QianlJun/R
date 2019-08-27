##p344
data(nutrient,package="flexclust")
head(nutrient,4)
d<-dist(nutrient)
as.matrix(d)[1:4,1:4]

##p346
data(nutrient,package="flexclust")
row.names(nutrient)<-tolower(row.names(nutrient))
nutrient.scaled<-scale(nutrient)
d<-dist(nutrient.scaled)
fit.average<-hclust(d,method="average")
plot(fit.average,hang=-1,cex=0.8,main="Average Linkage Clustering")

##p348
library(NbClust)
devAskNewPage(ask=TRUE)
nc<-NbClust(nutrient.scaled,distance="euclidean",
            min.nc=2,max.nc=15,method="average")
table(nc$Best.n[1,])
barplot(table(nc$Best.n[1,]),
        xlab="Numer of Clusters",ylab="Number of Criteria",
        main="Number of Clusters Chosen by 26 Criteria")

##获取最终的聚类方案
clusters<-cutree(fit.average,k=5)
table(clusters)
aggregate(nutrient,by=list(cluster=clusters),median)
aggregate(as.data.frame(nutrient.scaled),by=list(cluster=clusters),
          median)
plot(fit.average,hang=-1,cex=0.8,
     main="Average Linkage Clustering\n5 Cluster Solution")
rect.hclust(fit.average,k=5)

##351
wssplot<-function(data,nc=15,seed=1234){
  wss<-(nrow(data)-1)*sum(apply(data,2,var))
  for (i in 2:nc){
    set.seed(seed)
    wss[i]<-sum(kmeans(data,centers=i)$withiness)
  }
  plot(1:nc,wss,type="b",xlab="Number of Clusters",
       ylab="Within groups sum of squares")
}
data(wine,package="rattle")
head(wine)
library(rattle)
df<-scale(wine[-1])
wsplot