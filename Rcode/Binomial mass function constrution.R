choose(10,8)
choose(8,7)*0.5^7*0.5*1
(0.5^8)*8
choose(8,7)
Binomial
setwd("C:/Users/QJ/Desktop/workdata")
BFDeaths<-read.table("Birdfludeaths.txt",header = TRUE)
Deaths<-rowSums(BFDeaths[,2:16])
Deaths
names(Deaths)<-BFDeaths[,1]
Deaths
barplot(t(Counts),beside=TRUE)
Counts<-cbind(Cases,Deaths)
head(BFDeaths)
X1<-choose(5,1)*.85^1*.15^4
X1
X2<-choose(5,2)*0.85^2*.15^3
X3<-choose(5,3)*.85^3*.15^2
X4<-choose(5,4)*.85^4*.15^1
X5<-choose(5,5)*.85^5*.15^0
X5
Quantities<-c(X1,X2,X3,X4,X5
names(Quantities)<-c(1,2,3,4,5)
Quantities
barplot(Quantities,main="Probability of Different Quantities of Good Replicates")
setwd("C:/Users/QJ/Desktop/HW1/4_2")
