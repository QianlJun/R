data(sleep)
g1<-sleep$extra[1:10]
g2<-sleep$extra[11:20]
difference<-g2-g1
mn<-mean(difference) #1.67
s<-sd(difference) #1.13
n<-10
mn+c(-1,1)*qt(0.975,n-1)*s/sqrt(n)
t.test(difference)$conf.int
tStat<-sqrt(n)*mn/s
esVals<-seq(0,1,length=1000)
likVals<-dt(tStat,n-1,ncp=sqrt(n)*esVals)
likVals<-likVals/max(likVals)
plot(esVals,likVals,type="l")
lines(range(esVals[likVals>1/8]),c(1/8,1/8))
lines(range(esVals[likVals>1/16]),c(1/16,1/16))