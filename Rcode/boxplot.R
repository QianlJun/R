setwd("C:/Users/QJ/Desktop/HW1/4_2")
data1<-read.table("data_1.txt",header = TRUE)
data2<-read.table("data_2.txt",header = TRUE)
data1
mean(1,2,3)
mean(1,2,3)
a<-c(1,2,3)
mean(a)
mean(data1[],na,rm=TRUE)
data1
a
a<-data1
a
head(a)
mean(data1[2,2])
data1[,]
mean(data1[,])
data1mn<-mean(data1[,])
data1mn
data2mn<-mean(data2[,])
data2mn
boxplot(data1,main="Data1")
boxplot(data2,main="Data2")
histplot(data1,main="Data1")
?his
hist(data1[,])
?hist
barplot(data1[,])
hist(data1[,],prob=FALSE,main="Data1")
length(data1[,])
curve(dnorm,add=TRUE)
set.seed(123)
h<-rnorm(1000)
hist(h,prob=TRUE)
hist(data2[,],prob=TRUE,main="Data2")
data1_pvalue<-shapiro.test(data1[,])
data2_pvalue<-shapiro.test(data2[,])
qqplot(rnorm(length(data1[,])),data1[,],main="Q-Q Plot for Data1")
qqline(data1[,])
qqplot(rnorm(length(data2[,])),data2[,],main="Q-Q Plot for Data2")
qqline(data2[,])
factorial(6)
